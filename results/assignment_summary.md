# G3-F15 NOESY Assignment Summary

- Input XPK: `D:\structure_inference\noesy_g3f15x7.xpk`
- Sequence: `TCGGGAAGGGAGGG`
- Parsed peaks: 1245
- Best path score: 220.392
- Best-vs-second score gap: 0.211
- Policy: conservative; ambiguous states are reported instead of forced silent assignments.

## Best Assignment

| residue | base_ppm | h1prime_ppm | imino_ppm | intra_peak_ids | sequential_from_prev_peak_ids | confidence |
| --- | --- | --- | --- | --- | --- | --- |
| T1 | 7.5758 | 5.7824 |  | 42;214 |  | high |
| C2 | 7.7557 | 5.6465 |  | 47;49;235 | 41;233 | ambiguous |
| G3 | 7.6826 | 6.1940 | 11.1727 | 4;1017 | 51;1012;1225 | ambiguous |
| G4 | 7.7557 | 6.0739 | 11.1727 | 5;922 | 44;920 | high |
| G5 | 7.6826 | 5.9558 | 11.1727 | 43;1224 | 3;1015 | medium |
| A6 | 7.5758 | 5.8357 |  | 218 | 215 | medium |
| A7 | 8.0330 | 5.4221 |  | 0;560 | 7;559 | high |
| G8 | 7.8297 | 6.0739 | 11.1727 | 2;1045 | 1;1047 | ambiguous |
| G9 | 7.7557 | 5.0802 | 11.1727 | 872;924 | 5;922 | ambiguous |
| G10 | 7.9161 | 5.1320 | 11.1727 | 684;734 | 733;775 | ambiguous |
| A11 | 7.9737 | 6.1940 |  | 599;717 | 716 | ambiguous |
| G12 | 7.7557 | 5.4221 | 11.1727 | 48;232;254 | 44;920 | ambiguous |
| G13 | 7.8297 | 5.4640 | 11.1727 | 1046 | 1;1047 | ambiguous |
| G14 | 7.6826 | 5.5430 | 11.1727 | 1011 | 50;543 | high |

## Review Flags

| residue | confidence | notes |
| --- | --- | --- |
| C2 | ambiguous | top-path posterior 0.49 |
| G3 | ambiguous | top-path posterior 0.40 |
| G8 | ambiguous | H1prime cluster reused elsewhere in path; top-path posterior 0.63 |
| G9 | ambiguous | top-path posterior 0.14 |
| G10 | ambiguous | top-path posterior 0.14 |
| A11 | ambiguous | H1prime cluster reused elsewhere in path; top-path posterior 0.14 |
| G12 | ambiguous | H1prime cluster reused elsewhere in path; top-path posterior 0.49 |
| G13 | ambiguous | top-path posterior 0.32; weak intra/base-H1prime evidence |

## T Methyl Anchor Candidates

| low_ppm | high_ppm | peak_ids | intensity |
| --- | --- | --- | --- |
| 1.7523 | 7.4840 | 92 | 43.317 |
| 1.7523 | 7.5758 | 93;97 | 11.024 |
| 1.7523 | 7.3007 | 89;90;514 | 8.436 |
| 1.7523 | 7.6826 | 94;992 | 6.711 |
| 1.7523 | 7.7557 | 99 | 6.649 |
| 1.7523 | 7.6110 | 98 | 5.898 |
| 1.9095 | 7.3007 | 517 | 4.847 |
| 1.9095 | 7.6826 | 528 | 3.912 |
| 1.7523 | 7.8297 | 95 | 3.448 |
| 1.7523 | 7.9161 | 96;101 | 3.428 |

## Detected Shift Clusters

- H1prime/anomeric clusters: 5.011, 5.080, 5.132, 5.361, 5.422, 5.464, 5.543, 5.646, 5.693, 5.782, 5.836, 5.884, 5.956, 6.010, 6.074, 6.127, 6.194, 6.265, 6.317, 6.412
- Aromatic base clusters: 7.301, 7.484, 7.576, 7.611, 7.683, 7.756, 7.830, 7.916, 7.974, 8.033, 8.067
- Imino clusters: 10.901, 11.041, 11.173, 11.229, 12.066, 12.172, 12.308, 12.349, 12.416

## Rule Notes

- NOESY cross-peak intensity is used semi-quantitatively with the usual approximate r^-6 distance dependence; no absolute distances are claimed.
- H6/H8(i)-H1prime(i) and H6/H8(i)-H1prime(i-1) contacts are used for the sequential walk.
- T methyl-H6 support is used as the T1 anchor evidence.
- Guanine imino signals in the 10.5-12.5 ppm region are used as G-quadruplex validation evidence.

## References

- NOE intensity-distance approximation: https://pubs.rsc.org/en/content/articlehtml/2020/sc/d0sc02970j
- Automated NOE assignment and r^-6 volume relation: https://academic.oup.com/bioinformatics/article/23/3/381/235607
- G-quadruplex NOESY walking, thymine methyl, imino validation: https://academic.oup.com/nar/article/40/14/6946/2414847
