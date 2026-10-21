# Submission check — manuscript_v6.pdf and SupplementaryMaterial_v2.pdf

## A. Must fix before submission

1. **Table 1 (p. 7, l. 243) has lost every decimal point.** q, T_amb and T12−T8 read 453, 2978, −94 for 4.53 sL min⁻¹, 297.8 K, −9.4 K. Correct values (from `outputs/table_measured_envelope.md`):

| Run | G0 | q | T_amb | T̄w | T3 | T12−T8 |
|---|---:|---:|---:|---:|---:|---:|
| E81 | 256 | 4.53 | 297.8 | 700 | 550 | -9.4 |
| E80 | 256 | 6.61 | 296.6 | 686 | 561 | +25.8 |
| E79 | 256 | 8.04 | 296.7 | 669 | 560 | +39.6 |
| E78 | 256 | 10.01 | 298.1 | 643 | 554 | +54.0 |
| E77 | 256 | 13.85 | 297.1 | 589 | 525 | +58.1 |
| E76 | 304 | 4.53 | 296.4 | 908 | 647 | -97.9 |
| E75 | 304 | 6.95 | 296.7 | 911 | 681 | -29.1 |
| E74 | 304 | 9.03 | 296.5 | 896 | 697 | +13.8 |
| E73 | 304 | 13.17 | 296.9 | 840 | 694 | +56.6 |
| E72 | 304 | 18.32 | 296.4 | 757 | 653 | +60.7 |
| E71 | 456 | 7.13 | 300.8 | 1047 | 753 | -83.5 |
| E70 | 456 | 9.11 | 299.0 | 1045 | 779 | -31.1 |
| E69 | 456 | 10.50 | 296.7 | 999 | 770 | +11.5 |
| E68 | 456 | 12.50 | 298.7 | 979 | 773 | +35.1 |
| E67 | 456 | 15.28 | 296.3 | 935 | 764 | +53.5 |

2. **SM S6, closure sentence:** “1.148–1.323 at and 0.783–0.915 at 256” — the configuration “304” was dropped (the master has it). Should read “1.148–1.323 at 304 and 0.783–0.915 at 256”.
3. **Fend, Smirnova & Schöllgen (COMSOL contribution)** is cited as “(… Fend et al.)” with no year, and its reference entry has no year or venue. Add year 2011 and the venue (COMSOL Conference, DLR Institute of Solar Research) in EndNote so it prints as Fend et al., 2011 (or 2011a/b as EndNote assigns).
4. **SM S2 duplication:** the S2 paragraph repeats the sentence “The recurring output of … is a design recommendation — porosity, cell density, layering, channel aspect ratio …” twice in succession; and the sensor-selection result (121.6 vs 280.6 J K⁻¹, ratio 2.31) is stated twice in S3 (once in prose, once after the window table). Delete one of each.

## B. Grammar and typos (manuscript; line numbers as printed)

| Line | Text | Problem | Suggested fix |
|---:|---|---|---|
| 45 | He et al. (2023) states that the theoretical volumetric effect has not been validated by experiment | subject-verb disagreement: "et al." is treated as plural elsewhere in the text (e.g. "Capuano et al. (2016) illustrate", "Fend et al. (2013) reproduce", "Kribus et al. (2014) show") but singular here | change "states" to "state" |
| 108 | models agreeing with each other to within 25 K yet disagreeing with the experiment by 80–100 K signals | subject-verb disagreement: plural subject "models" paired with singular verb "signals" | change "signals" to "signal" |
| 113 | Second, presents an analysis framework (in section 3) that reduces the measurements | broken sentence: missing subject after "Second," | insert a subject, e.g. "Second, it presents an analysis framework" |
| 186 | Of which, three probes (type N, OD 1.5 mm, length 150 mm) contact the monolith side wall | Sentence fragment — "Of which," cannot begin an independent sentence | Change to "Of these, three probes..." or join to the preceding sentence with a comma |
| 203 | which is why we treat the reference-temperature choice as a systematic rather than a detail | Missing noun after "systematic" (e.g. "effect" or "uncertainty") makes the phrase ungrammatical | "as a systematic effect rather than a detail" |
| 212 | the entire expected monolith pressure drop, so the sensors readings are not useful | Wrong word form: "sensors readings" should be singular/possessive | "the sensor readings are not useful" |
| 214 | a corrective action has be devised since the quantity of pressure drop | "has be devised" is ungrammatical | "has to be devised" (or "has been devised") |
| 214 | can assist on discriminating between the two mechanisms discussed in section 5.2 | Wrong preposition: "assist on" should be "assist in" | "assist in discriminating" |
| 227 | they are retained in the dataset for repeatability check but excluded from every fit reported here | Missing article before "repeatability check" | "for a repeatability check" |
| 251 | Dependencies on additional, to the methodology framework assumptions, and other reference choices are declared in place. | Garbled sentence that cannot be clearly parsed | e.g. "Dependencies on additional methodology, framework assumptions, and other reference choices are declared in place." |
| 276 | appropriate for the heat transfer analysis herein but not explicit transport coefficients and transferable to another geometry | Missing negation makes the clause self-contradictory (should read "not... and not transferable") | "...but not explicit transport coefficients and not transferable to another geometry" |
| 290 | the ratio 𝑐𝑝/𝑘is being used rather calling it a Prandtl number | Missing word "than" after "rather" | "rather than calling it a Prandtl number" |
| 314 | So regressing the obtained λ against the gas advective conductance 𝑥 = 𝜀ṁ𝑐𝑝allow identification of both constants | Subject–verb disagreement: "regressing...allow" should be "allows" | "...allows identification of both constants" |
| 228 | All raw data are available in the work's repository (see Data Availability section) | Missing terminal period at end of sentence | "...(see Data Availability section)." |
| 412 | Where more than one reduction paths (referred as estimator choice) exists, as for the inversion | number disagreement: "one reduction paths" should be singular | "one reduction path ... exists" |
| 412 | reduction paths (referred as estimator choice) exists, as for the inversion threshold ε*, their | wrong/missing word: "referred as" should be "referred to as" | "referred to as estimator choice" |
| 511 | make a peripherally uniform flux the better idealization than a uniform temperature - in conjuction with the negligible web Biot number | misspelling: "conjuction" | "conjunction" |
| 421 | executed in the repository environment , with dependencies fixed by `Project.toml` and `Manifest.toml` | extra space before comma | "environment, with" |
| 364 | Three kinds of uncertainty are propagated and reported separately : i) Instrument error, ii) systematic bias | extra space before colon | "separately: i) Instrument" |
| 564 | a conductance fixed in position (i.e. only as function of length) has a fixed total over the monolith length | missing article: "as function" should be "as a function" | "(i.e. only as a function of length)" |
| 579 | This is one of main contributions of this work and further elaborated in section 5.2 | missing article: "one of main contributions" should be "one of the main contributions" | "This is one of the main contributions of this work" |
| 620 | Kribus et al., 2014; Moro Filho and Malalasekera, 2020) and of heterogeneous monolith models (Hayes et al., 2009 | wrong word: "built into X ... and of Y" breaks the parallel construction with "built into" | "...Malalasekera, 2020) and into heterogeneous monolith models (Hayes et al., 2009" |
| 834 | four models cluster within 25 K of each other and 80–100 K from experiment, mean-air-temperature measurement named a principal suspect | run-on/comma-spliced clauses with no connecting conjunction, cannot be parsed as a single sentence | split into separate sentences, e.g. "...80–100 K from experiment; mean-air-temperature measurement is named as a principal suspect..." |
| 857 | This is not an argument against reduced-order models but that the standard | missing word — should read "but an argument that" for parallel structure | This is not an argument against reduced-order models but an argument that the standard validation dataset cannot discriminate among them. |
| — | Their coefficient is volumetric and particle-referenced and ours apparent and wall-referenced | missing verb "is" before "apparent" | Their coefficient is volumetric and particle-referenced and ours is apparent and wall-referenced |
| 911 | two mechanisms have the potential to produce the required behaviour and sign, a transfer-unit count that rises with flow | garbled sentence — comma splice leaves "a transfer-unit count..." as an unparsable fragment | two mechanisms have the potential to produce the required behaviour and sign: a transfer-unit count that rises with flow. |
| 972 | rather than air enthalpy flow; second, the nearest published measurements isolate the effect | "second" has no corresponding "first" — broken parallel structure | insert "First," at the start of the preceding clause (e.g., "First, Broeske et al. (2023) state that...") |

## C. Checked and correct

- Every decimal number in the prose of both PDFs matches the pipeline archive or the markdown masters (0 unmatched across both files).
- Tables 2–5 and S3–S5: every numeric cell matches the generated output files.
- All 53 references in the list are cited in the text, and every author–year citation resolves to a list entry, including the a/b suffixes (Fend 2004a/b, Avila-Marin 2019a/b, Cornejo 2020a/b).
- Cornejo et al., Chem. Eng. Sci. 212 printing as 2020 is correct (volume 212 is the February 2020 issue); my earlier flag of 2020 vs 2019 is withdrawn.
- No leftover markup ($, \\, <mark>, **, EndNote temporary {…#…}, numeric [n], placeholders) in either PDF.
- SM cross-references to main-text sections (4.1, 4.2, 4.4, 4.7, 5.1–5.3) exist.
