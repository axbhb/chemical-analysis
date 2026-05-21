#!/usr/bin/env python3
"""NOESY-only assignment helper for G3-F15 peak lists.

The script reads an NMRView/CARA-style .xpk peak list, extracts evidence from
aromatic-anomeric, thymine-methyl, and guanine-imino regions, then performs a
conservative sequential-walk search over a supplied DNA sequence.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


AROMATIC_RANGE = (7.0, 8.6)
ANOMERIC_RANGE = (5.0, 6.5)
METHYL_RANGE = (1.45, 2.05)
IMINO_RANGE = (10.5, 12.5)

REFERENCES = [
    {
        "topic": "NOE intensity-distance approximation",
        "url": "https://pubs.rsc.org/en/content/articlehtml/2020/sc/d0sc02970j",
    },
    {
        "topic": "Automated NOE assignment and r^-6 volume relation",
        "url": "https://academic.oup.com/bioinformatics/article/23/3/381/235607",
    },
    {
        "topic": "G-quadruplex NOESY walking, thymine methyl, imino validation",
        "url": "https://academic.oup.com/nar/article/40/14/6946/2414847",
    },
]


@dataclass(frozen=True)
class Peak:
    idx: int
    hx: float
    hy: float
    hx_quality: str
    hy_quality: str
    volume: float
    intensity: float
    raw: str


@dataclass(frozen=True)
class OrientedPeak:
    idx: int
    low: float
    high: float
    volume: float
    intensity: float
    orientation: str


@dataclass
class Cluster:
    cid: int
    values: list[float]

    @property
    def mean(self) -> float:
        return sum(self.values) / len(self.values)

    @property
    def minimum(self) -> float:
        return min(self.values)

    @property
    def maximum(self) -> float:
        return max(self.values)

    def as_dict(self) -> dict[str, float | int]:
        return {
            "id": self.cid,
            "ppm": round(self.mean, 4),
            "n": len(self.values),
            "min": round(self.minimum, 4),
            "max": round(self.maximum, 4),
        }


@dataclass
class PairEvidence:
    low_cluster: int
    high_cluster: int
    low_ppm: float
    high_ppm: float
    peaks: list[OrientedPeak] = field(default_factory=list)

    @property
    def best_peak(self) -> OrientedPeak:
        return max(self.peaks, key=lambda peak: clean_number(peak.intensity))

    @property
    def max_intensity(self) -> float:
        return clean_number(self.best_peak.intensity)

    @property
    def total_intensity(self) -> float:
        return sum(clean_number(peak.intensity) for peak in self.peaks)

    @property
    def peak_ids(self) -> list[int]:
        return [peak.idx for peak in sorted(self.peaks, key=lambda peak: peak.idx)]


@dataclass(frozen=True)
class CandidateState:
    residue_type: str
    sugar_cluster: int
    base_cluster: int
    sugar_ppm: float
    base_ppm: float
    intra: PairEvidence
    type_score: float
    methyl: PairEvidence | None = None
    imino: PairEvidence | None = None
    imino_shift: float | None = None
    notes: tuple[str, ...] = ()

    @property
    def state_key(self) -> str:
        return f"s{self.sugar_cluster}:b{self.base_cluster}"


@dataclass
class WalkPath:
    score: float
    states: list[CandidateState]
    edges: list[PairEvidence | None]
    notes: list[list[str]]


def clean_number(value: float) -> float:
    if value is None or math.isnan(value) or math.isinf(value):
        return 0.0
    return float(value)


def in_range(value: float, bounds: tuple[float, float]) -> bool:
    return bounds[0] <= value <= bounds[1]


def parse_float(token: str) -> float:
    try:
        return float(token)
    except ValueError:
        return float("nan")


def parse_xpk(path: Path) -> list[Peak]:
    peaks: list[Peak] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if len(parts) < 18 or not parts[0].isdigit():
            continue
        peaks.append(
            Peak(
                idx=int(parts[0]),
                hx=parse_float(parts[2]),
                hy=parse_float(parts[9]),
                hx_quality=parts[5],
                hy_quality=parts[12],
                volume=parse_float(parts[15]),
                intensity=parse_float(parts[16]),
                raw=line,
            )
        )
    if not peaks:
        raise ValueError(f"No peak rows could be parsed from {path}")
    return peaks


def orient_peak(
    peak: Peak,
    low_range: tuple[float, float],
    high_range: tuple[float, float],
) -> OrientedPeak | None:
    if in_range(peak.hx, low_range) and in_range(peak.hy, high_range):
        return OrientedPeak(
            peak.idx, peak.hx, peak.hy, peak.volume, peak.intensity, "Hx-low/Hy-high"
        )
    if in_range(peak.hy, low_range) and in_range(peak.hx, high_range):
        return OrientedPeak(
            peak.idx, peak.hy, peak.hx, peak.volume, peak.intensity, "Hy-low/Hx-high"
        )
    return None


def collect_oriented(
    peaks: Iterable[Peak],
    low_range: tuple[float, float],
    high_range: tuple[float, float],
) -> list[OrientedPeak]:
    oriented: list[OrientedPeak] = []
    for peak in peaks:
        found = orient_peak(peak, low_range, high_range)
        if found is not None:
            oriented.append(found)
    return oriented


def make_clusters(values: Iterable[float], tol: float) -> list[Cluster]:
    clusters: list[Cluster] = []
    for value in sorted(clean_number(v) for v in values if not math.isnan(v)):
        if not clusters or abs(value - clusters[-1].mean) > tol:
            clusters.append(Cluster(len(clusters), [value]))
        else:
            clusters[-1].values.append(value)
    return clusters


def cluster_id(value: float, clusters: list[Cluster], tol: float) -> int | None:
    best_id: int | None = None
    best_delta = float("inf")
    for cluster in clusters:
        delta = abs(value - cluster.mean)
        if delta < best_delta:
            best_id = cluster.cid
            best_delta = delta
    if best_id is None or best_delta > tol:
        return None
    return best_id


def build_pair_map(
    oriented: Iterable[OrientedPeak],
    low_clusters: list[Cluster],
    high_clusters: list[Cluster],
    tol: float,
) -> dict[tuple[int, int], PairEvidence]:
    pair_map: dict[tuple[int, int], PairEvidence] = {}
    for peak in oriented:
        low_id = cluster_id(peak.low, low_clusters, tol)
        high_id = cluster_id(peak.high, high_clusters, tol)
        if low_id is None or high_id is None:
            continue
        key = (low_id, high_id)
        if key not in pair_map:
            pair_map[key] = PairEvidence(
                low_cluster=low_id,
                high_cluster=high_id,
                low_ppm=low_clusters[low_id].mean,
                high_ppm=high_clusters[high_id].mean,
            )
        pair_map[key].peaks.append(peak)
    return pair_map


def ppm_window_score(value: float, preferred: tuple[float, float], broad: tuple[float, float]) -> float:
    if in_range(value, preferred):
        return 2.0
    if in_range(value, broad):
        return 0.5
    nearest = min(abs(value - broad[0]), abs(value - broad[1]))
    return -1.0 - min(nearest * 3.0, 3.0)


def residue_type_score(residue_type: str, sugar_ppm: float, base_ppm: float) -> tuple[float, list[str]]:
    notes: list[str] = []
    if residue_type == "T":
        score = ppm_window_score(base_ppm, (7.25, 7.85), (7.05, 8.05))
        score += ppm_window_score(sugar_ppm, (5.45, 5.90), (5.20, 6.00))
    elif residue_type == "C":
        score = ppm_window_score(base_ppm, (7.25, 7.90), (7.00, 8.10))
        score += ppm_window_score(sugar_ppm, (5.20, 5.90), (5.00, 6.20))
    elif residue_type == "G":
        score = ppm_window_score(base_ppm, (7.45, 8.20), (7.20, 8.35))
        score += ppm_window_score(sugar_ppm, (5.20, 6.30), (5.00, 6.45))
    elif residue_type == "A":
        score = ppm_window_score(base_ppm, (7.30, 8.35), (7.00, 8.55))
        score += ppm_window_score(sugar_ppm, (5.20, 6.35), (5.00, 6.50))
    else:
        raise ValueError(f"Unsupported residue type: {residue_type}")

    if score < 1.0:
        notes.append("outside preferred shift window")
    return score, notes


def find_best_high_cluster_support(
    support_map: dict[int, list[PairEvidence]],
    high_cluster: int,
) -> PairEvidence | None:
    items = support_map.get(high_cluster, [])
    if not items:
        return None
    return max(items, key=lambda evidence: evidence.max_intensity)


def build_support_by_high_cluster(
    pair_map: dict[tuple[int, int], PairEvidence],
) -> dict[int, list[PairEvidence]]:
    grouped: dict[int, list[PairEvidence]] = {}
    for (_low, high), evidence in pair_map.items():
        grouped.setdefault(high, []).append(evidence)
    return grouped


def make_imino_support(
    imino_aromatic: list[OrientedPeak],
    imino_clusters: list[Cluster],
    base_clusters: list[Cluster],
    tol: float,
) -> dict[int, list[PairEvidence]]:
    pair_map = build_pair_map(imino_aromatic, imino_clusters, base_clusters, tol)
    return build_support_by_high_cluster(pair_map)


def candidate_states(
    sequence: str,
    pair_map: dict[tuple[int, int], PairEvidence],
    methyl_support: dict[int, list[PairEvidence]],
    imino_support: dict[int, list[PairEvidence]],
    min_intensity: float,
    max_candidates_per_type: int,
) -> dict[str, list[CandidateState]]:
    candidates: dict[str, list[CandidateState]] = {base: [] for base in sorted(set(sequence))}
    for (_sugar, _base), evidence in pair_map.items():
        if evidence.max_intensity < min_intensity:
            continue
        for residue_type in candidates:
            score, notes = residue_type_score(residue_type, evidence.low_ppm, evidence.high_ppm)
            score += math.log1p(evidence.max_intensity)

            methyl = None
            if residue_type == "T":
                methyl = find_best_high_cluster_support(methyl_support, evidence.high_cluster)
                if methyl is None:
                    score -= 3.0
                    notes.append("no thymine methyl support")
                else:
                    score += 2.5 + math.log1p(methyl.max_intensity)

            imino = None
            imino_shift = None
            if residue_type == "G":
                imino = find_best_high_cluster_support(imino_support, evidence.high_cluster)
                if imino is not None:
                    score += 1.5 + min(math.log1p(imino.max_intensity), 3.5)
                    imino_shift = imino.low_ppm
                else:
                    notes.append("no imino-aromatic support")

            candidates[residue_type].append(
                CandidateState(
                    residue_type=residue_type,
                    sugar_cluster=evidence.low_cluster,
                    base_cluster=evidence.high_cluster,
                    sugar_ppm=evidence.low_ppm,
                    base_ppm=evidence.high_ppm,
                    intra=evidence,
                    type_score=score,
                    methyl=methyl,
                    imino=imino,
                    imino_shift=imino_shift,
                    notes=tuple(notes),
                )
            )

    for residue_type, items in candidates.items():
        items.sort(key=lambda item: item.type_score, reverse=True)
        candidates[residue_type] = items[:max_candidates_per_type]
    return candidates


def edge_score(
    previous: CandidateState,
    current: CandidateState,
    pair_map: dict[tuple[int, int], PairEvidence],
) -> tuple[float, PairEvidence | None, list[str]]:
    notes: list[str] = []
    edge = pair_map.get((previous.sugar_cluster, current.base_cluster))
    if edge is None:
        notes.append("missing sequential H6/H8(i)-H1prime(i-1) peak")
        return -8.0, None, notes

    score = 2.0 + 1.8 * math.log1p(edge.max_intensity)
    if edge.max_intensity < 2.5:
        score -= 1.5
        notes.append("weak sequential peak")
    if previous.sugar_cluster == current.sugar_cluster:
        score -= 6.0
        notes.append("same H1prime cluster as previous residue")
    if previous.base_cluster == current.base_cluster:
        if previous.residue_type != current.residue_type:
            score -= 7.0
            notes.append("same aromatic cluster as previous different residue type")
        else:
            score -= 2.0
            notes.append("same aromatic cluster as previous residue")
    if edge.best_peak.idx == previous.intra.best_peak.idx:
        score -= 5.0
        notes.append("sequential evidence reuses previous intra peak")
    if edge.best_peak.idx == current.intra.best_peak.idx:
        score -= 3.0
        notes.append("sequential evidence reuses current intra peak")
    return score, edge, notes


def search_walk(
    sequence: str,
    candidates: dict[str, list[CandidateState]],
    pair_map: dict[tuple[int, int], PairEvidence],
    beam_size: int,
) -> list[WalkPath]:
    beams: list[WalkPath] = [WalkPath(score=0.0, states=[], edges=[], notes=[])]
    for residue_index, residue_type in enumerate(sequence):
        next_beams: list[WalkPath] = []
        residue_candidates = candidates.get(residue_type, [])
        if not residue_candidates:
            raise ValueError(f"No candidates found for residue type {residue_type}")

        for path in beams:
            used_states = {state.state_key for state in path.states}
            for candidate in residue_candidates:
                state_penalty = 0.0
                candidate_notes = list(candidate.notes)
                if candidate.state_key in used_states:
                    state_penalty -= 18.0
                    candidate_notes.append("exact state reused in path")
                sugar_reuse_count = sum(
                    1 for state in path.states if state.sugar_cluster == candidate.sugar_cluster
                )
                base_reuse_count = sum(
                    1 for state in path.states if state.base_cluster == candidate.base_cluster
                )
                if sugar_reuse_count:
                    state_penalty -= 2.25 * sugar_reuse_count
                    candidate_notes.append("H1prime cluster reused elsewhere in path")
                if base_reuse_count > 1:
                    state_penalty -= 0.75 * (base_reuse_count - 1)
                if path.states:
                    edge_delta, edge, edge_notes = edge_score(path.states[-1], candidate, pair_map)
                    score = path.score + candidate.type_score + edge_delta + state_penalty
                    next_beams.append(
                        WalkPath(
                            score=score,
                            states=path.states + [candidate],
                            edges=path.edges + [edge],
                            notes=path.notes + [candidate_notes + edge_notes],
                        )
                    )
                else:
                    start_bonus = 2.0 if candidate.methyl is not None else -2.0
                    score = path.score + candidate.type_score + start_bonus + state_penalty
                    next_beams.append(
                        WalkPath(
                            score=score,
                            states=[candidate],
                            edges=[],
                            notes=[candidate_notes],
                        )
                    )

        next_beams.sort(key=lambda path: (round(path.score, 6), path_preference(path, sequence)), reverse=True)
        beams = next_beams[:beam_size]
        if not beams:
            raise ValueError(f"Beam search failed at residue {residue_index + 1}{residue_type}")
    return beams


def path_preference(path: WalkPath, sequence: str) -> float:
    """Tie-break only: favor chemically plausible starts without hiding ambiguity."""
    preference = 0.0
    if not path.states:
        return preference
    first = path.states[0]
    if sequence[0] == "T":
        if first.methyl is not None:
            preference += 2.0
        if 5.50 <= first.sugar_ppm <= 5.82:
            preference += 1.5
        if 7.25 <= first.base_ppm <= 7.80:
            preference += 1.0
    if len(path.states) > 1 and sequence[1] == "C":
        second = path.states[1]
        if 7.25 <= second.base_ppm <= 7.90:
            preference += 2.0
        if 5.20 <= second.sugar_ppm <= 5.75:
            preference += 1.0
    return preference


def intensity_class(intensity: float) -> str:
    if intensity >= 20.0:
        return "strong"
    if intensity >= 6.0:
        return "medium"
    if intensity >= 2.0:
        return "weak"
    return "very_weak"


def weighted_state_frequencies(paths: list[WalkPath], top_n: int) -> list[dict[str, float]]:
    selected = paths[:top_n]
    if not selected:
        return []
    best = selected[0].score
    weights = [math.exp(max(path.score - best, -50.0)) for path in selected]
    total = sum(weights)
    freqs: list[dict[str, float]] = []
    for residue_index in range(len(selected[0].states)):
        counts: dict[str, float] = {}
        for path, weight in zip(selected, weights):
            key = path.states[residue_index].state_key
            counts[key] = counts.get(key, 0.0) + weight / total
        freqs.append(counts)
    return freqs


def confidence_for(
    state: CandidateState,
    edge: PairEvidence | None,
    residue_notes: list[str],
    frequency: float,
) -> tuple[str, list[str]]:
    notes = list(residue_notes)
    confidence = "high"
    if frequency < 0.65:
        confidence = "ambiguous"
        notes.append(f"top-path posterior {frequency:.2f}")
    elif frequency < 0.85:
        confidence = "medium"
        notes.append(f"top-path posterior {frequency:.2f}")

    if state.intra.max_intensity < 2.5:
        confidence = "low_confidence" if confidence != "ambiguous" else confidence
        notes.append("weak intra/base-H1prime evidence")
    elif state.intra.max_intensity < 6.0 and confidence == "high":
        confidence = "medium"
        notes.append("weak intra/base-H1prime evidence")
    if edge is None and state.residue_type != "T":
        confidence = "low_confidence" if confidence == "high" else confidence
    elif edge is not None and edge.max_intensity < 2.5:
        confidence = "low_confidence" if confidence == "high" else confidence
    elif edge is not None and edge.max_intensity < 6.0 and confidence == "high":
        confidence = "medium"
        notes.append("weak sequential evidence")
    if state.residue_type == "G" and state.imino is None:
        notes.append("G imino not directly linked to assigned aromatic cluster")
    return confidence, notes


def evidence_ids(evidence: PairEvidence | None) -> str:
    if evidence is None:
        return ""
    return ";".join(str(idx) for idx in evidence.peak_ids[:8])


def evidence_intensity(evidence: PairEvidence | None) -> str:
    if evidence is None:
        return ""
    return f"{evidence.max_intensity:.3f}"


def path_to_rows(best_path: WalkPath, all_paths: list[WalkPath], sequence: str, top_n: int) -> list[dict[str, str]]:
    freqs = weighted_state_frequencies(all_paths, top_n)
    rows: list[dict[str, str]] = []
    for index, state in enumerate(best_path.states):
        edge = best_path.edges[index - 1] if index else None
        freq = freqs[index].get(state.state_key, 0.0) if freqs else 1.0
        confidence, notes = confidence_for(state, edge, best_path.notes[index], freq)
        rows.append(
            {
                "residue_index": str(index + 1),
                "residue": f"{sequence[index]}{index + 1}",
                "residue_type": sequence[index],
                "base_ppm": f"{state.base_ppm:.4f}",
                "h1prime_ppm": f"{state.sugar_ppm:.4f}",
                "imino_ppm": "" if state.imino_shift is None else f"{state.imino_shift:.4f}",
                "intra_peak_ids": evidence_ids(state.intra),
                "intra_intensity": evidence_intensity(state.intra),
                "intra_class": intensity_class(state.intra.max_intensity),
                "sequential_from_prev_peak_ids": evidence_ids(edge),
                "sequential_from_prev_intensity": evidence_intensity(edge),
                "methyl_peak_ids": evidence_ids(state.methyl),
                "methyl_intensity": evidence_intensity(state.methyl),
                "imino_peak_ids": evidence_ids(state.imino),
                "imino_intensity": evidence_intensity(state.imino),
                "state_score": f"{state.type_score:.3f}",
                "path_state_frequency": f"{freq:.3f}",
                "confidence": confidence,
                "notes": "; ".join(dict.fromkeys(notes)),
            }
        )
    return rows


def summarize_path(path: WalkPath, sequence: str) -> dict[str, object]:
    residues = []
    for index, state in enumerate(path.states):
        edge = path.edges[index - 1] if index else None
        residues.append(
            {
                "residue": f"{sequence[index]}{index + 1}",
                "base_ppm": round(state.base_ppm, 4),
                "h1prime_ppm": round(state.sugar_ppm, 4),
                "imino_ppm": None if state.imino_shift is None else round(state.imino_shift, 4),
                "intra_peak_ids": state.intra.peak_ids,
                "sequential_from_prev_peak_ids": [] if edge is None else edge.peak_ids,
                "methyl_peak_ids": [] if state.methyl is None else state.methyl.peak_ids,
                "imino_peak_ids": [] if state.imino is None else state.imino.peak_ids,
                "notes": path.notes[index],
            }
        )
    return {"score": round(path.score, 4), "residues": residues}


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join("---" for _ in columns) + " |"
    body = []
    for row in rows:
        body.append("| " + " | ".join(row.get(column, "") for column in columns) + " |")
    return "\n".join([header, separator] + body)


def top_support_rows(support: dict[int, list[PairEvidence]], limit: int = 10) -> list[dict[str, str]]:
    evidences = [item for items in support.values() for item in items]
    evidences.sort(key=lambda evidence: evidence.max_intensity, reverse=True)
    rows = []
    for evidence in evidences[:limit]:
        rows.append(
            {
                "low_ppm": f"{evidence.low_ppm:.4f}",
                "high_ppm": f"{evidence.high_ppm:.4f}",
                "peak_ids": evidence_ids(evidence),
                "intensity": evidence_intensity(evidence),
            }
        )
    return rows


def write_summary(
    path: Path,
    args: argparse.Namespace,
    peaks: list[Peak],
    best_path: WalkPath,
    rows: list[dict[str, str]],
    all_paths: list[WalkPath],
    sugar_clusters: list[Cluster],
    base_clusters: list[Cluster],
    imino_clusters: list[Cluster],
    methyl_support: dict[int, list[PairEvidence]],
) -> None:
    ambiguous = [row for row in rows if row["confidence"] in {"ambiguous", "low_confidence"}]
    score_gap = ""
    if len(all_paths) > 1:
        score_gap = f"{all_paths[0].score - all_paths[1].score:.3f}"

    lines = [
        "# G3-F15 NOESY Assignment Summary",
        "",
        f"- Input XPK: `{Path(args.xpk).resolve()}`",
        f"- Sequence: `{args.sequence}`",
        f"- Parsed peaks: {len(peaks)}",
        f"- Best path score: {best_path.score:.3f}",
        f"- Best-vs-second score gap: {score_gap or 'n/a'}",
        f"- Policy: conservative; ambiguous states are reported instead of forced silent assignments.",
        "",
        "## Best Assignment",
        "",
        markdown_table(
            rows,
            [
                "residue",
                "base_ppm",
                "h1prime_ppm",
                "imino_ppm",
                "intra_peak_ids",
                "sequential_from_prev_peak_ids",
                "confidence",
            ],
        ),
        "",
        "## Review Flags",
        "",
    ]
    if ambiguous:
        lines.append(
            markdown_table(
                ambiguous,
                ["residue", "confidence", "notes"],
            )
        )
    else:
        lines.append("No low-confidence or ambiguous residues were flagged by the current thresholds.")

    lines.extend(
        [
            "",
            "## T Methyl Anchor Candidates",
            "",
            markdown_table(top_support_rows(methyl_support), ["low_ppm", "high_ppm", "peak_ids", "intensity"]),
            "",
            "## Detected Shift Clusters",
            "",
            f"- H1prime/anomeric clusters: {', '.join(f'{c.mean:.3f}' for c in sugar_clusters)}",
            f"- Aromatic base clusters: {', '.join(f'{c.mean:.3f}' for c in base_clusters)}",
            f"- Imino clusters: {', '.join(f'{c.mean:.3f}' for c in imino_clusters)}",
            "",
            "## Rule Notes",
            "",
            "- NOESY cross-peak intensity is used semi-quantitatively with the usual approximate r^-6 distance dependence; no absolute distances are claimed.",
            "- H6/H8(i)-H1prime(i) and H6/H8(i)-H1prime(i-1) contacts are used for the sequential walk.",
            "- T methyl-H6 support is used as the T1 anchor evidence.",
            "- Guanine imino signals in the 10.5-12.5 ppm region are used as G-quadruplex validation evidence.",
            "",
            "## References",
            "",
        ]
    )
    for reference in REFERENCES:
        lines.append(f"- {reference['topic']}: {reference['url']}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_json(
    path: Path,
    args: argparse.Namespace,
    peaks: list[Peak],
    sugar_clusters: list[Cluster],
    base_clusters: list[Cluster],
    imino_clusters: list[Cluster],
    paths: list[WalkPath],
    sequence: str,
    top_n: int,
) -> None:
    payload = {
        "metadata": {
            "input_xpk": str(Path(args.xpk).resolve()),
            "sequence": sequence,
            "parsed_peak_count": len(peaks),
            "ppm_tolerance": args.ppm_tol,
            "min_intensity": args.min_intensity,
            "policy": "conservative",
            "references": REFERENCES,
        },
        "clusters": {
            "h1prime_anomeric": [cluster.as_dict() for cluster in sugar_clusters],
            "aromatic_base": [cluster.as_dict() for cluster in base_clusters],
            "imino": [cluster.as_dict() for cluster in imino_clusters],
        },
        "top_paths": [summarize_path(path, sequence) for path in paths[:top_n]],
    }
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    xpk_path = Path(args.xpk)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sequence = args.sequence.upper().replace("U", "T")
    unsupported = sorted(set(sequence) - set("ACGT"))
    if unsupported:
        raise ValueError(f"Unsupported residue letters in sequence: {unsupported}")

    peaks = parse_xpk(xpk_path)
    aromatic_anomeric = collect_oriented(peaks, ANOMERIC_RANGE, AROMATIC_RANGE)
    methyl_base = collect_oriented(peaks, METHYL_RANGE, AROMATIC_RANGE)
    imino_aromatic = collect_oriented(peaks, IMINO_RANGE, AROMATIC_RANGE)

    sugar_clusters = make_clusters((peak.low for peak in aromatic_anomeric), args.ppm_tol)
    base_clusters = make_clusters((peak.high for peak in aromatic_anomeric), args.ppm_tol)
    imino_clusters = make_clusters((peak.low for peak in imino_aromatic), args.ppm_tol)

    aromatic_pair_map = build_pair_map(
        aromatic_anomeric, sugar_clusters, base_clusters, args.ppm_tol * 1.5
    )
    methyl_pair_map = build_pair_map(methyl_base, make_clusters((p.low for p in methyl_base), args.ppm_tol), base_clusters, args.ppm_tol * 1.5)
    methyl_support = build_support_by_high_cluster(methyl_pair_map)
    imino_support = make_imino_support(
        imino_aromatic, imino_clusters, base_clusters, args.ppm_tol * 1.5
    )

    candidates = candidate_states(
        sequence,
        aromatic_pair_map,
        methyl_support,
        imino_support,
        args.min_intensity,
        args.max_candidates_per_type,
    )
    paths = search_walk(sequence, candidates, aromatic_pair_map, args.beam_size)
    paths.sort(key=lambda path: (round(path.score, 6), path_preference(path, sequence)), reverse=True)
    best_path = paths[0]
    rows = path_to_rows(best_path, paths, sequence, args.top_n)

    csv_path = out_dir / "assignment_table.csv"
    json_path = out_dir / "assignment_candidates.json"
    summary_path = out_dir / "assignment_summary.md"

    write_csv(csv_path, rows)
    write_json(
        json_path,
        args,
        peaks,
        sugar_clusters,
        base_clusters,
        imino_clusters,
        paths,
        sequence,
        args.top_n,
    )
    write_summary(
        summary_path,
        args,
        peaks,
        best_path,
        rows,
        paths,
        sugar_clusters,
        base_clusters,
        imino_clusters,
        methyl_support,
    )

    print(f"Parsed peaks: {len(peaks)}")
    print(f"Aromatic/anomeric peaks: {len(aromatic_anomeric)}")
    print(f"T methyl/base peaks: {len(methyl_base)}")
    print(f"Imino/aromatic peaks: {len(imino_aromatic)}")
    print(f"Best path score: {best_path.score:.3f}")
    print(f"Wrote: {csv_path}")
    print(f"Wrote: {json_path}")
    print(f"Wrote: {summary_path}")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Conservative NOESY-only assignment helper for DNA G-quadruplex peak lists."
    )
    parser.add_argument("--xpk", default="noesy_g3f15x7.xpk", help="Input .xpk peak list")
    parser.add_argument("--sequence", default="TCGGGAAGGGAGGG", help="DNA sequence")
    parser.add_argument("--out-dir", default="results", help="Output directory")
    parser.add_argument("--ppm-tol", type=float, default=0.03, help="Chemical-shift clustering tolerance")
    parser.add_argument("--min-intensity", type=float, default=1.5, help="Minimum peak intensity for candidate states")
    parser.add_argument("--top-n", type=int, default=10, help="Number of candidate paths to retain in JSON/report scoring")
    parser.add_argument("--beam-size", type=int, default=5000, help="Beam size for sequential-walk search")
    parser.add_argument("--max-candidates-per-type", type=int, default=80, help="Candidate states retained for each residue type")
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
