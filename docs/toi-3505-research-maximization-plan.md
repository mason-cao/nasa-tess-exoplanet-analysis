# TOI-3505.01 maximum-research blueprint

| Field | Value |
|---|---|
| Working project | Ground-based follow-up and dilution-aware multi-epoch vetting of TOI-3505.01 |
| Working paper title | **Stress-Testing TOI-3505.01: A Dilution-Aware, Multi-Epoch Analysis Using TESS, Gaia, and Public TFOP Constraints** (the final title must match the achieved result) |
| Target | TOI-3505.01 / TIC 390988385 |
| Document status | Living research design, not a statement that every analysis has already been completed |
| Last updated | 2026-07-20 (personal-strategy, repository, and AERIS-methods audit integrated; external target facts remain on the dated 2026-07-14 snapshot unless noted) |
| Program fit | GMU NASA Data Science & Astronomy Research Internship, Summer 2026 |
| Hard deadlines | Symposium poster/talk: Aug 1, 2026 · Paper due: Aug 29, 2026 |
| Strategic identity | Evidence engineering for noisy physical systems: domain-correct astronomy first, reproducible research software second |

## What the pivot re-research found

This document replaces the TOI-6241.01 blueprint after the project pivot to TOI-3505.01. Every external fact below was re-verified against live archives (NASA Exoplanet Archive TAP, ExoFOP-TESS JSON, Gaia DR3 TAP, MAST/TESScut, exo.MAST TCE index) on 2026-07-14, and the delivered GMU FITS files were opened and their headers read. The findings are not cosmetic — several of them overturn assumptions embedded in the pivot request itself, and the plan below is built around them. In order of severity:

1. **The delivered GMU night does not contain a transit under the current ephemeris.** The 283 delivered science frames span BJD_TDB 2459782.598–2459782.839 (2022-07-22 02:13–08:00 UT, 5.77 h). Propagating the current catalog ephemeris (P = 2.9151556 ± 0.0000117 d, T0 = 2459793.534385 ± 0.0021) backwards, the nearest transits fall at BJD ≈ 2459781.8738 and ≈ 2459784.7889. The night covers orbital phases **0.249–0.331** — no primary transit, no secondary eclipse (phase 0.5) either. The miss is ~0.72 d against a back-propagated timing uncertainty of ~3 minutes; this is not a marginal call. The likeliest explanation is that the night was scheduled against a stale 2021 ephemeris (the TOI was created 2021-06-23, presumably from Sector 14 data taken in 2019, and its parameters were not updated to the current values until 2024-08-22), but that is a hypothesis to test, not a fact. The GMU analysis is therefore reframed: identity + photometric quality + eclipse *upper limit at the covered phases* + ephemeris archaeology, not transit recovery. Do not soften this. Lecture 8 explicitly anticipates nights that do not show the transit — see [Lecture 8, result possibilities](../data_and_lectures/Lecture8_Schar2025_clean_redesign.pdf#page=14).
2. **A fourth TESS sector exists.** The pivot thesis names Sectors 14, 41, and 54. MAST/TESScut additionally lists **Sector 81** (2024, camera 2, CCD 3), and exo.MAST lists TCEs for `s0054`, `s0081`, and a **multi-sector run spanning s0014–s0086**. A "multi-epoch consistency" thesis that silently ignores the most recent epoch and the existing pipeline Data Validation products would be indefensible under review. The thesis below is extended to all four sectors (2019–2024), with S14/41/54 retained as the named core.
3. **This is a genuinely crowded, genuinely blended target — the "dilution-aware" framing is not decoration, it is the whole problem.** The field sits at galactic latitude −3.46°. The TIC v8.2 contamination ratio is **0.547** (contaminating flux ≈ 55% of target flux in the TESS aperture model). Worse: SOAR speckle imaging (2021-10-01) shows a companion at **0.517″ with ΔI = 1.7 mag** (~21% of the target's flux) that has **no separate Gaia DR3 entry** — meaning it is absent from the TIC and therefore absent from *every* catalog-driven dilution correction QLP or SPOC applied. Shane/ShARCS AO (2021-07-19) catalogs 11 companions between 0.507″ and 8.355″. No seeing-limited photometry — TESS, GMU, or any public SG1 light curve — resolves the 0.52″ pair.
4. **The star itself is not settled.** Gaia DR3 gives the target (source 1824561377891719424) **RUWE = 3.69**, far beyond the ~1.4 threshold where the single-star astrometric model has failed — unsurprising given the 0.52″ companion, but it poisons the parallax (DR3: 1.245 mas, naively ~800 pc) which conflicts with the TIC distance (373.5 pc, DR2-era). TRES reconnaissance spectroscopy on ExoFOP reports **Teff 6423 ± 63 K, log g 3.875 ± 0.109, v sin i 14.74 km/s** — a possibly evolved, fast-rotating F star — against TIC's 6220 K / 1.335 R☉ dwarf-ish solution. The candidate radius (7.46 R⊕ per the TOI table) is a hostage of these unresolved choices: defensible stellar-radius scenarios move it between roughly 7.5 and ~18 R⊕ (screening arithmetic in the stellar section). The candidate is either a rare hot-Neptune-desert object, an ordinary hot Jupiter, or an eclipsing binary — **and dilution plus stellar characterization decide which**. That is the paper.
5. **The public TFOP record is rich — which kills "first follow-up" novelty and replaces it with something better.** ExoFOP lists **7 time-series observations** (KeplerCam ip full; ULMT rp full *with an NEB check*; MuSCAT2 g,r,i,z_s egress — four simultaneous colors; CMO g′ full; OAA Ic full; and a **prior GMU 0.8 m R-band night, 2021-06-28, tagged by the mentor**), 2 imaging campaigns (Shane AO J/Ks, SOAR speckle I), 3 TRES spectra, and **Keck/HIRES: "6 RVs between UT 2022-09-01 and 2022-09-12"**. ~147 files are downloadable with a free ExoFOP login. Nobody has published any of it: the NASA Archive disposition is still PC, `pscomppars` has no entry for this TIC, and exact-identifier searches return no dedicated paper (novelty ledger entries #1–4). The defensible contribution is therefore the *synthesis and stress-test*: a student-led, dilution-aware, multi-epoch reanalysis that unifies four TESS sectors, the public SG1 archive, Gaia, and the imaging constraints — plus one GMU night (ours) that appears nowhere on ExoFOP.
6. **The GMU night is simultaneous with TESS Sector 54.** S54 ran ~2022-07-09 to 2022-08-05 (verify against the data-release notes); the GMU night of 2022-07-21/22 sits inside it. Same-night space photometry of the same star is a free, rare cross-check on the GMU night's photometric fidelity — and it converts an off-transit night into a calibrated precision benchmark.
7. **Program-calendar integration is unchanged.** Same Summer 2026 program, same lecture rules, same deadlines (Discord light-curve gate → symposium Aug 1 → paper Aug 29). All compliance sections are carried forward intact and re-pointed at this target.

## Purpose

This document is the complete research blueprint for making the TOI-3505.01 project as rigorous, original, and scientifically useful as the available data allow while staying inside the internship rules.

The required paper remains a ground-based TESS Object of Interest follow-up study: determine what the GMU observation of TOI-3505.01 does and does not show at the predicted time, depth, and duration — which, for this night, includes determining *why the predicted transit was not in the observing window* and what the night constrains anyway. Every advanced component below strengthens that central result; none replaces it.

The project can produce a defensible paper from exactly the situation we are in: a nondetection-by-geometry night, an ephemeris audit, and a multi-epoch space-based consistency analysis. The program explicitly expects imperfect nights: many are "plagued by clouds, data gaps" and end inconclusive, marginal, partial, or false-positive, and "only a few will result in clear, on target detections" — see [Lecture 8, result possibilities](../data_and_lectures/Lecture8_Schar2025_clean_redesign.pdf#page=14). What the program does not tolerate — and what this document is designed to prevent — is claims that outrun the evidence.

## Student strategy and project identity — INTERNAL, not paper text

The personal value of this project is not that exoplanets can be relabeled as environmental science. The authentic connection is methodological: both this study and Project AERIS ask how to make reliable claims from noisy, incomplete, differently sampled observations of a physical system. The astronomy must remain domain-native; the cross-project throughline is **evidence verification**, not "AI applied to two unrelated topics."

A dated audit of the public AERIS repository on 2026-07-20 found a mature verification architecture: frozen evaluation sets linked to code and data hashes, measurement-process channel grouping, explicit non-independence caveats, leave-one-source/channel-out ablations, trigger-channel circularity tests, blinded expert-label packets, provenance gates, numeric claim checks, and a large automated test suite. This TOI project should transfer those research habits where they genuinely apply:

- Freeze analysis decisions before seeing the answer.
- Distinguish a new measurement from a new reduction of the same measurement.
- Track the provenance and dependence of every evidence stream.
- Preserve supporting, contradicting, silent, and unusable evidence rather than forcing a binary verdict.
- Use ablations, negative controls, synthetic recovery, and human review to learn what actually carries a conclusion.
- Set a claim ceiling and abstain when the available data cannot separate scenarios.

The transfer stops there. AERIS is not evidence about TOI-3505.01, its terminology does not belong in the astronomy paper, and no LLM or multi-agent layer should be added to make the projects look connected. The strongest portfolio story is that the same student can learn the physics and measurement process of a new domain, then build trustworthy software around it.

### Current public-evidence gap (2026-07-20)

The NASA repository currently proves an AU Mic/TESS introductory analysis, a substantial TOI-3718 practice reduction, and a TOI-3505 archive/header/timing audit with progress figures. It does **not yet** prove a completed TOI-3505 differential light curve, four-sector reanalysis, Bayesian or hierarchical inference, journal manuscript, or automated test suite. The public description must stay at "analyzing" or "building" until those artifacts exist. Closing this evidence gap is more important than strengthening the wording.

Conversely, AERIS already contains extensive evaluation and testing infrastructure but its public README still labels the formal freeze, expert labels, and final statistics as unfinished. Do not use this internship to imply that AERIS has a completed result it does not yet have, or use AERIS's code volume as evidence of TOI-3505 progress. Each project earns its own claims.

### Two legitimate academic readings

| Audience | What this project can honestly demonstrate | What not to claim |
|---|---|---|
| Earth systems / environmental engineering | Remote-sensing literacy; reconciliation of space- and ground-based measurements; calibration, spatial resolution, contamination, temporal sampling, uncertainty, and physical-model constraints | That exoplanet research is environmental research, or that astronomical dilution directly models atmospheric source attribution |
| Computer science | Reproducible scientific pipelines; configuration-driven experiments; tests for time conversion and dilution algebra; provenance and data lineage; statistically valid evaluation; clear human-review gates | That software sophistication substitutes for an astronomical result, or that using more models makes the work more scientific |
| Shared throughline | Reliable inference from heterogeneous measurements of real physical systems | A generic "AI + NASA + climate" brand |

### Decision rule for maximizing the internship

Add a task only if it does at least one of the following without endangering the paper deadline:

1. changes or tightens the scientific claim ceiling;
2. tests whether a result depends on a fragile analysis choice;
3. creates a reusable, mentor-verifiable research artifact;
4. teaches a domain method that can be explained without scripted language.

If a task mainly makes the project sound more advanced, postpone it. Admissions-facing summaries, portfolio copy, and public case studies are downstream products created after the result is frozen; they never decide which analysis is run or which result is emphasized.

## Scope legend

Each proposed analysis is assigned one of five levels.

| Level | Meaning | Rule |
|---|---|---|
| **CORE** | Needed for the internship paper | Complete before pursuing advanced claims. |
| **HIGH-VALUE BONUS** | Feasible with public or existing data and likely to improve the paper | Prioritize after the first trustworthy light curve. |
| **ADVANCED** | Statistically or technically demanding | Use only if the method, assumptions, and code can be explained and validated. |
| **CONDITIONAL** | Useful only for certain outcomes or if additional data exist | Activate only when its dependency is satisfied. |
| **GATED / FUTURE** | Requires mentor approval, new observations, restricted products, or a separate paper | Describe as future work unless approval and data are obtained. |

## One-sentence research strategy

Determine whether our analysis can localize the 2.9151556-day signal attributed to TOI-3505.01 to the unresolved TIC 390988385 system (or show that it belongs to a resolvable neighbor), and whether its depth, duration, shape, and timing remain physically consistent across TESS Sectors 14, 41, 54, and 81 once dilution, aperture selection, stellar variability, and false-positive scenarios are explicitly modeled — using the GMU 2022-07-22 night as an identity-verified photometric benchmark and ephemeris-audit anchor, and the public TFOP archive (multi-band SG1 light curves, NEB checks, AO/speckle imaging, TRES spectroscopy) as quantitative constraints rather than as an appeal to authority.

## Strongest defensible novelty thesis

The pivot request states the thesis as:

> Can TOI-3505.01's 2.915-day signal be independently localized to the target star and shown to remain physically consistent across TESS Sectors 14, 41, and 54 after accounting for dilution, aperture selection, stellar variability, and false-positive scenarios?

This is the right *kind* of thesis for this target — but as stated it has two defects that a hostile reviewer would find in five minutes, so the executed version must repair them:

1. **It omits Sector 81 and the existing SPOC/TESS-SPOC Data Validation products.** Four sectors spanning 2019–2024 exist, along with TCEs (and therefore pipeline centroid/difference-image/odd-even diagnostics) for s0054, s0081, and a multi-sector s0014–s0086 run. "Independently localized" claims that ignore the pipeline's own localization tests are not independent — they are incomplete. The executed thesis covers all four sectors and treats DV reproduction as part of the localization evidence.
2. **"Localized to the target star" may be unachievable as worded, and the plan must say so up front.** The 0.517″ ΔI = 1.7 companion is unresolved by TESS (21″ pixels), by GMU (~1–3″ seeing), and by every public SG1 light curve. Centroid and NEB analyses can localize the signal to *the unresolved pair*, not to the primary alone. Separating the pair requires chromatic-depth arguments (the public multi-band SG1 data make this genuinely testable), transit-derived stellar-density consistency, AO-resolved photometry, or RVs. The claim ladder below caps what each evidence layer can support.

The executed thesis: **a student-led, dilution-aware stress test of TOI-3505.01 across all four TESS epochs (2019–2024) and the public TFOP archive, quantifying (a) where the signal originates at each achievable angular scale, (b) whether one physically self-consistent transit solution survives contact with all the public evidence simultaneously, and (c) what the candidate is if it does not.**

The paper should aim to contribute five evidence layers:

1. **Localization evidence:** per-sector difference images, centroid offsets, and aperture-dependence of depth, cross-checked against pipeline DV products and the public SG1 NEB checks, with explicit statements of the angular scales actually excluded — and the 0.52″ pair explicitly *not* excluded unless chromatic/density evidence earns it.
2. **Dilution evidence:** a decomposed contamination budget (cataloged neighbors via TIC/TESS-cont, plus the uncataloged speckle companion), applied consistently to every depth from every instrument, with the candidate radius reported under each defensible stellar scenario rather than as one falsely precise number.
3. **Multi-epoch consistency evidence:** a predeclared comparison of per-sector (and, where comparable, per-instrument) depths, durations, and shapes over 2019–2024, testing whether one physical transit model explains all epochs or whether the between-epoch scatter demands systematics, variability, or a false-positive explanation. A hierarchical Bayesian model is the preferred advanced implementation only after the simpler per-sector estimates, dependence map, and injection yardstick are valid.
4. **Stellar and variability evidence:** rotation/pulsation characterization from the TESS time series (v sin i predicts P_rot/sin i ≲ 5 d — squarely measurable in 27-d sectors), spectroscopic-vs-transit-derived stellar density, and an honest treatment of the RUWE-poisoned distance.
5. **Timing evidence:** a refined linear ephemeris across 2019–2024, an O−C diagram, the reconstruction of why the 2022 GMU night missed the transit, and concrete 2026 re-observation windows.

No document or abstract should call the work "first," "novel," "confirmed," or "validated" until the corresponding checks in this blueprint are completed and the mentor approves the wording. The working title therefore uses "Analysis," not "Validation"; see the result-dependent title section.

## Program-scope compliance and deliverables calendar

### What the program requires

The lecture materials frame the paper around GMU campus-telescope follow-up of a TESS candidate and explicitly allow inconclusive results. The central test is whether a transit appears on the correct target near the expected time with a compatible depth and duration, using a finder chart to confirm target identity. For this night, the honest execution of that test is: verify identity, reduce the night to a trustworthy light curve, show quantitatively that the predicted transit (current ephemeris) was not in the window, determine what ephemeris *was* used at scheduling, and report what the light curve constrains. See [Lecture 8, paper science goals](../data_and_lectures/Lecture8_Schar2025_clean_redesign.pdf#page=13), [Lecture 8, result possibilities](../data_and_lectures/Lecture8_Schar2025_clean_redesign.pdf#page=14), and [the paper-structure guide](../data_and_lectures/Paperstructure.pdf).

### Fixed program dates (Summer 2026)

| Date | Deliverable | Notes |
|---|---|---|
| Ongoing | Discord progress posts | Post light curves, AIJ settings, measurement tables (`.xls`), and `.plotcfg` files so results are reproducible by others; see [Lecture 8, paper logistics](../data_and_lectures/Lecture8_Schar2025_clean_redesign.pdf#page=10). |
| After first light curve posted | TFOP-context unlock | Only Dr. Plavchan and Kevin can access internal TESS-team information (timing uncertainty, TFOP notes); they share it *after* a light curve is posted. For this target the TFOP notes plausibly include the HIRES RV outcome and SG1 dispositions — posting early is a scientific dependency of the first order, not an administrative step. |
| July 18, 25 | Live office hours | Bring reduction questions, the ephemeris-archaeology question (what ephemeris was used to schedule 2022-07-21?), and NEB/dilution questions here. |
| July 31 | Symposium tours | — |
| **Aug 1** | **Symposium: poster and/or talk** | Poster 48×36″ max (40×30″ template provided); talk 5 min max, PPT/Google Slides, slides follow the paper sections. A survey selects poster vs talk; talk selection depends on analysis posted to Discord. |
| **Aug 29** | **Paper due (PDF)** | Send to execed@gmu.edu and nasa.schar.program@gmail.com. Word, Google Docs, or Overleaf/LaTeX (aastex). Papers are lightly edited and posted as a GMU student journal. |
| Aug 25 | Mentor availability drops | GMU semester starts; support tapers. Front-load anything needing mentor input — especially TFOP-note access and claim-ceiling review. |

### Collaboration and authorship rules from the lectures

- Teams of up to 3 may analyze the same data set; the sign-up spreadsheet records this. **Audit the spreadsheet for the TOI-3505.01 row(s):** who signed up, what predicted ingress/egress times were listed (this is primary evidence for the ephemeris-archaeology question), and whether anyone else is analyzing this night.
- The data-collecting observers should be offered co-authorship so the paper can say "we" for the observations. The FITS headers credit the "Omegalambda automation code" as observer and Dr. Plavchan as software owner — identify the responsible humans for the 2022-07-21 night via the spreadsheet/mentor and credit them; Dr. Plavchan is an author with the GMU Department of Physics and Astronomy affiliation. See [Paperstructure, Observations](../data_and_lectures/Paperstructure.pdf#page=11).
- Swap drafts and compare reduction approaches on Discord; self-organized study groups are encouraged.
- Do not pay for access to any paper; request copies through the mentor or use arXiv preprints.

### What the lectures permit as advanced work

The paper-structure slides describe statistical false-positive work such as `vespa` as beyond the normal internship scope but welcome as a bonus, with the stated ultimate goal of ground-based follow-up being a <1% false-positive probability built from ExoFOP imaging, spectroscopy, and light-curve tests. That makes advanced validation permissible, not mandatory — and it must not be attempted at the expense of a correct reduction. See [the advanced-analysis note](../data_and_lectures/Paperstructure.pdf#page=16). For this target, note the sobering corollary: the same slide's <1% goal may be *unreachable* while the 0.52″ companion scenario stands untested by any resolved observation — plan for a claim ceiling below "validated" and treat anything better as upside.

The lectures also explicitly name the AstroImageJ **nearby eclipsing binary (NEB) check** as a program-supported tool the mentors will help run and interpret — see [Lecture 8, NEB note](../data_and_lectures/Lecture8_Schar2025_clean_redesign.pdf#page=13). Two public NEB checks (ULMT 2021-10-15, KeplerCam 2022-06-14) already exist for this target; the program-aligned move is to run our own on our own night *and* re-derive theirs from their posted measurement tables, stating clearly that an out-of-transit night constrains neighbor variability, not the transit's source.

### What is gated

JWST work is a later opportunity requiring approval after the campus-telescope light curve is approved, and it is a separate paper on a separate observation. It appears only in the gated extension section. See [Lecture 8, JWST gate](../data_and_lectures/Lecture8_Schar2025_clean_redesign.pdf#page=4).

### AI and authorship rule

Any code, equation, model, or wording used in the final project must be understood and defensible by the student: "Do not submit text, claims, or analysis you do not understand. Use AI cautiously and transparently." Preserve an analysis log and disclose assistance. See [Lecture 1, AI guidance](../data_and_lectures/Lecture1_intro.pdf#page=26). This applies with special force to the Bayesian machinery below: a hierarchical model you cannot explain at office hours is a liability, not a bonus.

For every AI-assisted contribution retained in the project, the student must be able to explain the input data, algorithm, assumptions, failure modes, and at least one test that could falsify the output. Keep a compact assistance ledger (`date`, `tool`, `task`, `student verification`, `retained/changed/rejected`). Do not publish generated wording as a scientific claim until it has been checked against the claim–evidence ledger described below.

## Target fact sheet

Values below were re-verified against the live NASA Exoplanet Archive TOI table, ExoFOP-TESS JSON, Gaia DR3 TAP, MAST/TESScut, and exo.MAST on **2026-07-14**. The TOI parameter set was last updated 2024-08-22 (source `qlp-s54-tois`); the TOI itself was created 2021-06-23 ("found in faint-star QLP search"). Re-snapshot immediately before the final paper.

| Quantity | Value (2026-07-14) | Interpretation or caution |
|---|---:|---|
| TIC ID | 390988385 | Stable target identifier to use in every search. |
| TOI | 3505.01 | Candidate identifier. |
| Candidate provenance | QLP faint-star search; TOI created 2021-06-23; parameters from `qlp-s54-tois` (updated 2024-08-22) | The 2021 creation date vs 2024 parameter update is the seed of the ephemeris-archaeology question. |
| TESS / TFOP disposition | PC / PC | Planet candidate; not confirmed or validated by this label — despite 5 years of TFOP activity, which is itself informative. |
| Right ascension | 19:48:10.43 / 297.043476° | Matches the plate-solve targets and the delivered FITS headers (RAOBJ2K 19.802897 h). |
| Declination | +18:41:56.1 / +18.698914° | Matches FITS headers (DECOBJ2K +18.698914°). |
| Galactic latitude | −3.46° | Galactic plane. Crowding is the defining property of this field: >1000 Gaia DR3 sources within 2.5′. |
| TESS magnitude | 10.9372 ± 0.0067 | TIC v8.2. Bright enough for strong ground-based SNR; check saturation at 50 s on the GMU 0.8 m. |
| V magnitude | 11.243 ± 0.011 | — |
| Gaia G | 11.476 (DR3, this query) vs 11.3355 ± 0.003 (ExoFOP listing) | The discrepancy is itself evidence of blended, DR-dependent photometry; record both with provenance. |
| J / H / K | 10.346 / 10.101 / 10.056 | 2MASS; blended (0.52″ pair unresolved). |
| Period | 2.9151556 ± 0.0000117 d | TOI ephemeris; starting point, not fixed truth. |
| Reference epoch | BJD 2459793.534385 ± 0.0020787 | Confirm time system is BJD_TDB in the source; epoch sits inside Sector 54 (2022-08-02). |
| Duration | 2.004 ± 0.21 h | Short for P = 2.915 d around an R ≥ 1.3 R☉ star — implies high impact parameter under every stellar scenario (screening estimate b ≈ 0.78 for R* = 1.34 R☉; ≥0.9 if subgiant). A load-bearing red flag: grazing geometries are the classic EB disguise. Refit rather than inherit. |
| Depth | 2910 ± 196 ppm | 3.164 ± 0.213 mmag. Determine what dilution correction QLP applied before treating this as a physical depth. |
| Radius ratio implied by depth | ≈ 0.0539 | `sqrt(0.00291)`, before any dilution correction. |
| Candidate radius | 7.4578 R⊕ (no error published) | Hostage to stellar radius and dilution; see the scenario table in the stellar section. Do not quote it bare. |
| SNR (pipeline) | 19 | — |
| Host Teff | 6220 K (TIC) vs **6423 ± 63 K (TRES)** | Late-F. |
| Host radius | 1.335 R☉ (TIC) | Derived from a DR2-era distance; suspect (see RUWE row). |
| Host log g | **3.875 ± 0.109 (TRES)**; TIC value absent | Subgiant-suggestive. If real, the candidate radius grows by ~60%+. |
| Host [m/H] | 0.183–0.298 ± 0.08 (three TRES epochs) | Metal-rich. |
| v sin i | 14.74 ± 0.5 km/s (TRES) | Fast rotator → P_rot/sin i ≲ 4.6 d (R*=1.34 R☉) to ~7.6 d (2.2 R☉). Rotation signal should be *visible in TESS* — a direct time-series test of the stellar scenarios. |
| Distance | 373.477 pc (TIC, DR2-era) vs Gaia DR3 parallax 1.2448 mas (naively ~800 pc) | **Both suspect**: RUWE 3.69 invalidates the single-star astrometric solution. Distance is a dominant systematic on R*, hence on Rp. |
| Insolation / predicted mass / predicted RV | 504.7 S⊕ / 43.44 M⊕ / 17.2 m/s | Model predictions, not measurements; the predicted RV assumes the Neptune-ish radius scenario. |
| TSM / ESM | 86.9 / 21.5 | Prioritization metrics, not detections. |
| TIC contamination ratio | **0.547163** | ~55% of the target's flux worth of *cataloged* contaminating flux in the TESS aperture model — and the 0.52″ companion is NOT in this number (no Gaia/TIC entry). |
| Speckle companion | **0.5171″, ΔI = 1.7 (SOAR HRCam, 2021-10-01)**; contrast Δ6.6 mag @ 1″ | Flux ratio ≈ 0.21. Unresolved by TESS, GMU, and all seeing-limited SG1 photometry. The central false-positive threat. |
| AO companions | 11 sources, 0.507″–8.355″, Δ(J/Ks) 1.53–7.13 (Shane/ShARCS, 2021-07-19, PI Dressing) | Download the annotated PDFs and sensitivity curves; these feed TRICERATOPS directly. |
| Gaia DR3 source | 1824561377891719424 | G = 11.476, BP−RP = 0.745, **RUWE = 3.694**, parallax 1.2448 mas, pmRA +4.72 / pmDec +2.95 mas/yr, `non_single_star = 0`, no variability flag. RUWE ≫ 1.4 is astrometric evidence that the single-star solution is poor and is compatible with unresolved multiplicity; note `non_single_star = 0` merely means no NSS solution was fit, not that the star is single. ExoFOP's PM (2 ± 1.3, 1.7 ± 1.3) disagrees with DR3 — record both, trust neither blindly. |
| Nearest Gaia neighbors | 5.70″ (G 20.05), 5.73″ (G 20.53), 6.33″ (G 19.21), 6.58″ (G 18.93), 6.88″ (G 18.21), 8.42″ (G 17.52) … | All bright-enough-to-matter screening arithmetic is in the blend section; inside 30″ only the G 17.52 star could formally mimic the diluted signal (requiring a ~76% eclipse). |
| TESS coverage | **Sectors 14, 41, 54, 81** (TESScut: cam 1/1/2/2, CCD 3/3/4/3) | ExoFOP still lists "14,41,54" — S81 verified via TESScut on 2026-07-14. FFI cadences: 30 min (S14), 10 min (S41, S54), 200 s (S81). Whether 2-min target data exist in any sector must be verified on MAST (lightkurve search) — do not assume either way. |
| TCEs (exo.MAST) | `s0014-s0086: TCE_1`, `s0054-s0054: TCE_1`, `s0081-s0081: TCE_1` | Pipeline Data Validation products exist — retrieve the DV reports/mini-reports/time series for all three and identify which pipeline (SPOC 2-min vs TESS-SPOC FFI) produced each. |
| Public follow-up (ExoFOP) | **7 time series, 2 imaging campaigns, 5 spectroscopy entries, ~147 files, 0 CTOIs** | Specifics in the next section. Counts alone are not constraints; download and inspect the products. |
| Future TESS coverage | Unknown — WTV was erroring on 2026-07-14 | Re-check the Web TESS Viewing tool / tess-point for Cycle 8+ coverage; a 2026 northern sector would add a live epoch. |

Primary catalog starting points: [NASA Exoplanet Archive TOI table](https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=TOI) and [ExoFOP-TESS target page](https://exofop.ipac.caltech.edu/tess/target.php?id=390988385).

## Public follow-up assets (verified 2026-07-14)

Each row is a concrete product to obtain, read, and use. ExoFOP files are downloadable with a free ExoFOP login; confirm citation/acknowledgment rules with the mentor before using them in the paper (SG1 products are community contributions — credit observers by name, and honor any usage notes such as OAA's "Please consult notes file").

| Asset | Facility / instrument | Date | Key content | Use in this project |
|---|---|---|---|---|
| Time series (prior GMU night) | **GMU 0.8 m / SBIG STX-16803, R** | 2021-06-28 | "Full" coverage, 272 points, 2 runs, tag: plavchan; 11 files | An in-house prior epoch, tagged full-transit. Reanalyze from the posted measurement tables; ask the mentor whether the raw FITS exist on the program drives — a two-epoch GMU comparison (one on-transit, one off-transit) is a free major upgrade. |
| Time series + NEB check | ULMT 0.6 m / STX-16803, r′ | 2021-10-15 | "Transit and NEB check"; 14 files incl. NEB results | Independent localization evidence at seeing-limited scales; re-derive the clearance verdicts from the posted tables rather than citing the tag. |
| Time series + NEB check | FLWO 1.2 m / KeplerCam, i′ | 2022-06-14 | Full coverage, 741 points; 10 files incl. NEBcheck tar | Highest-point-count public light curve; second independent NEB check. |
| Time series (multicolor) | TCS 1.52 m / **MuSCAT2, g,r,i,z_s simultaneous** | 2023-07-14 | Egress coverage, 429 points; 34 files incl. transit fits | **The chromaticity test.** Blend/EB scenarios are generically chromatic; four simultaneous bands constrain the depth-vs-wavelength slope even with partial coverage. |
| Time series | CMO SAI RC-600 0.6 m, g′ | 2023-05-05 | Full coverage, 143 points | Blue-end depth point for the chromatic analysis. |
| Time series | OAA 0.4 m, Ic | 2023-07-11 | Full coverage, 100 points; "Please consult notes file" | Red-end depth point; read the notes file before use. |
| AO imaging | Shane 3 m / ShARCS, J + Ks | 2021-07-19 | 11 companions 0.507″–8.355″, Δmag 1.53–7.13 (PI Dressing) | Companion census + contrast curves; direct TRICERATOPS input; two-band Δmag gives a crude color (hence temperature/size guess) for the 0.52″ companion. |
| Speckle imaging | SOAR 4.1 m / HRCam, I | 2021-10-01 | Companion at 0.5171″, ΔI = 1.7; sensitivity Δ6.6 @ 1″ | Independent confirmation of the close companion in an optical band; sensitivity curve for validation input. |
| Reconnaissance spectra ×3 | FLWO 1.5 m / TRES, R = 44,000 | 2021-08-02, 2021-10-04, 2021-10-07 | SNR ≈ 29–32; derived Teff/log g/[m/H]/v sin i on ExoFOP; classification plots + extracted FITS | Three epochs can reveal large RV shifts (SB1) and composite spectra (SB2); the log g = 3.875 value drives the subgiant scenario — scrutinize its uncertainty and methodology before leaning on it. |
| Precision RVs (existence known, data restricted) | Keck 10 m / HIRES | 2022-08-08 (1 epoch); **2022-09-01→12 (6 RVs)** | Listed on ExoFOP; products not public | Someone measured an orbit-scale RV dataset in 2022 and has not published. Ask the mentor for TFOP-note context *after* posting the light curve. Also a novelty risk: a paper could appear mid-project. |
| Time-series photometry of the 2022-07-21 GMU night | — | — | **Not on ExoFOP** | Our night is publicly undocumented — the one dataset this project alone brings to the table. |

## The ephemeris problem — CORE, and the night's most interesting result

### The verified geometry

- Delivered window (mid-exposure BJD_TDB, from headers; re-derive independently before publication): **2459782.59823 → 2459782.83867** (5.77 h, 283 × 50 s frames, ~73.7 s effective cadence, no large gaps).
- Current-ephemeris transits nearest the window: **T0 − 4P = 2459781.87376** (2022-07-21 08:58 UT) and **T0 − 3P = 2459784.78892** (2022-07-24 06:56 UT).
- Phase coverage of the window: **0.2485 → 0.3310**. Half-duration in phase is ~0.014; the window is ~0.24 in phase away from any transit and ~0.17 from the phase-0.5 secondary. Back-propagated 1σ timing uncertainty at the window epoch: ~3.0 minutes. The miss is unambiguous *under the current ephemeris*.

### The archaeology tasks

1. Reconstruct the TOI ephemeris as it stood in mid-2022 (the alert-era values): the sign-up spreadsheet's listed period and predicted ingress/egress for this night are primary evidence; the mentor and TFOP notes may retain the alert parameters; the Swarthmore transit-finder inputs used at scheduling may be recoverable.
2. Independently refit Sector 14 alone (the only data existing at TOI creation in June 2021) and propagate that ephemeris, with covariance, to 2022-07-22. Quantify: does an S14-only ephemeris predict a transit inside the GMU window? If yes, the night is a clean, quantified case study of stale-ephemeris drift — a known, citable failure mode of TOI follow-up (see the ephemeris-refinement literature in the references) and a genuinely useful program-level result.
3. If the S14-only reconstruction does *not* put a transit in the window, escalate honestly: candidate explanations include scheduling from a different tool/ephemeris, a deliberate out-of-transit baseline request, or plain error. Ask at office hours; document the answer.
4. Deliver both ephemerides side by side: the prediction as it plausibly stood at scheduling, and the current best prediction. The difference *is* a publishable figure.

### What the off-transit night is still good for

- **Identity and photometric benchmark:** plate-solved, calibrated, BJD-audited photometry of the field, with precision quantified against the simultaneous Sector 54 data.
- **Eclipse upper limit at phases 0.25–0.33:** a planet-on-target scenario predicts a flat light curve; a detected ppt-level dip would be evidence *against* the clean interpretation (off-ephemeris events, alias-period EBs, variability). Fit for it under the frozen protocol; report the limit.
- **Neighbor variability monitoring:** light curves of every NEB-check star that night (variables, flares, deep eclipses at other phases are all findings).
- **Stellar variability:** 5.8 h of R-band monitoring of a v sin i ≈ 15 km/s F star, simultaneous with TESS — a cross-instrument variability check.
- **Ephemeris-drift case study:** item 2 above.

Do not let anyone — including yourself — describe this night as a failed observation. It is an off-transit observation with a reconstructable cause and several quantitative uses. But equally: do not describe it as a transit observation. It is not one, under the current ephemeris, and the paper must say so in plain sentences.

## Intended GMU observation (as delivered)

| Field | Value (from FITS headers and file inventory; spreadsheet row still to be audited) |
|---|---|
| Night | 2022-07-21/22 UT (local evening of July 21) |
| Science frames | 283 × 50.000 s, filter "Red" (R), numbered 0001–0283, no missing numbers |
| Time span | 02:13:25 → 07:59:37 UT (start times); mid-exposure BJD_TDB 2459782.59823 → 2459782.83867 |
| Airmass run | 1.32 → 1.066 (culmination between frames ~140–210) → 1.40 |
| Calibration delivered | 10 × R flats (3.5 s, filenames tagged `-final` — determine what processing produced them before trusting), 10 × 50 s darks, 10 × 3.5 s darks (flat-matched), 11 unique focus FITS + 22 focus PNGs. **No bias frames** — matched-exposure darks make bias separable only if the pipeline is documented; record the adopted calibration algebra explicitly. |
| Telescope | GMU 0.8 m (APTDIA 812.8 mm), f = 5769 mm |
| Camera | SBIG STX-16803, 4096², 9 μm pixels, gain 1.31 e−/ADU, −20.4 °C setpoint (stable at −20.37 in headers) |
| Plate scale / FOV | 0.322″/px, 22.0′ × 22.0′ |
| Site | 38:49:41.4 N, 77:18:19.2 W, 154 m (GMU campus observatory) |
| Header time products | JD (start), JD_UTC (≈ mid-exposure), BJD_TDB (mid-exposure), HJD — treat all as *claims to verify*, not truth; re-derive BJD_TDB with astropy from DATE-OBS + site + coordinates and reconcile |
| Acquisition | "Omegalambda automation code" (OBSERVER header); MaxIm DL 6.17; FLIPSTAT 'Mirror/Rotate 90 CW' |
| Meridian passage | Between frames ~140 and ~210 — inspect for pointing/rotation/PSF discontinuities and treat as a candidate systematics breakpoint in detrending |

Planning checks (not observed counts): 283 frames × 50 s = 3.93 h open shutter over a 5.77 h span (68% shutter duty cycle; ~24 s/frame overhead — consistent with full-frame readout of a 16.8 Mpix CCD). Saturation check required: Tmag 10.94 / V 11.24 at 50 s on a 0.8 m in good seeing may approach full well in the core — measure peak ADU on early frames before choosing apertures.

Also audit the program's 2021–2025 data-drive folders and the full sign-up spreadsheet for any *additional* GMU nights on TOI-3505: **at least one exists** (2021-06-28, on ExoFOP, tagged full coverage). Recovering its raw frames in-house would give a two-epoch GMU dataset — one on-transit, one off-transit — which is a materially stronger paper.

## Data-identity protocol — CORE

Target identity is the first scientific validation test of any night: the program itself requires confirming the transit occurs "on the TESS object of interest" using a finder chart. The header RA/Dec already match the TIC position to the arcsecond, but headers are claims, not proof — plate-solve. Before any interpretation, complete this standard protocol on the delivered files:

Build a manifest with one row per file:

- Original path and filename; SHA-256 checksum (compute at extraction; the six delivered zip archives total ~6.3 GB — decide deliberately whether extracted FITS live inside the repo, and if so behind Git LFS with a quota check, or outside it with a manifest; do not casually commit 6 GB).
- FITS object/target name; DATE-OBS; exposure time; filter; image type; binning; detector temperature; observatory coordinates when present.
- Header right ascension and declination.
- Plate-solved field center and angular separation from TOI-3505.01.
- Science, flat, dark, bias, test, focus, or unknown classification.
- Keep/reject decision and a written reason.
- Any filename-versus-header disagreement.

Identity acceptance requires the available evidence to agree: observing log, object headers, plate solution, finder chart ([Swarthmore finder-chart tool](https://astro.swarthmore.edu/transits/finding_charts.cgi)), target position in the field, date, filter, and exposure sequence. A filename alone is insufficient; a spreadsheet row alone is insufficient. In this galactic-plane field (>1000 Gaia sources in 2.5′), misidentifying the target among near-equal neighbors is a live risk — overlay the Gaia catalog on the plate solution and confirm the target by coordinates, not by appearance. Record the identity verdict in the analysis ledger before reducing.

## Novelty verification protocol

**Status as of 2026-07-14** (novelty ledger entries #1–5):

1. Exact-name searches for "TOI-3505," "TOI 3505," "TOI-3505 b," and "TIC 390988385" reveal **no dedicated publication**.
2. The NASA Exoplanet Archive planetary-systems table (`pscomppars`) has **no entry** for this TIC — the candidate is not a confirmed planet under any name.
3. ExoFOP lists **no CTOIs** and extensive unpublished TFOP products, including Keck/HIRES RVs from 2022 — a strong hint that a professional team (the HIRES time suggests a hot-Jupiter/giant-candidate RV program) holds unpublished results. **This is the single largest novelty risk: a confirmation/FP paper could appear at any time.** Set a weekly ADS/arXiv alert on the TIC and TOI identifiers.
4. The GMU Schar Astro-Scholars journal archive (JASR, three volumes) contains **no TOI-3505 paper**. Caution: a **TOI-3506.01** paper exists — adjacent number; do not let searches, Discord threads, or spreadsheet lookups conflate them.
5. The RAVEN FFI survey (sectors 1–55, P = 0.5–16 d) plausibly *searched* this star — S14/41/54 are all ≤ 55 and the period is in range — but membership and classification are **unverified**. Check the released tables (Zenodo record 19661443) for TIC 390988385 before writing anything about RAVEN and this target.

Consequences for wording:

- **Not safe:** "first ground-based follow-up" (seven public SG1 time series exist), "no prior analysis" (two NEB checks and pipeline DV reports exist), "validated" (nothing supports it yet).
- **Safe now:** "TOI-3505.01 has remained a planet candidate since 2021 despite an extensive unpublished TFOP record; we present a dilution-aware synthesis of public evidence, a student-led multi-epoch analysis of all four TESS sectors, and a previously undocumented GMU observation." Recheck "previously undocumented" against program drives and public archives before submission.
- **Potential after a clean re-audit and mentor approval:** sharper versions, if the audit holds.

Before final drafting, repeat the audit:

1. Search NASA ADS by TOI number, TIC ID, coordinates, and host aliases (Gaia DR3 1824561377891719424, 2MASS/UCAC4 aliases from ExoFOP).
2. Search arXiv with the same identifiers and spelling variants ("TOI 3505", "TOI-3505b").
3. Check the current NASA Exoplanet Archive disposition and reference fields, and `pscomppars` again.
4. Check current ExoFOP time-series, imaging, spectroscopy, stellar, notes tabs, and any new CTOIs.
5. Ask the mentor whether unpublished TFOP products, student projects, or papers in preparation cover this target — specifically the HIRES program's status.
6. Check the JASR archive again for a TOI-3505 paper.
7. Record every query, exact search string, date, result URL, and export in the novelty ledger; repeat once more before submission.

## Research questions

### Primary localization question — CORE

Can the 2.9151556-d signal be shown, from per-sector TESS difference images, centroid offsets, aperture-dependence of depth, the pipeline DV products, and the public SG1 NEB checks, to originate within the unresolved TIC 390988385 system (target + 0.52″ companion) rather than from any resolvable neighbor — with quantified angular exclusion limits at each step?

### Dilution question — CORE

What is the dilution-corrected transit depth in each sector and each public instrument under a single declared contamination model (TIC/TESS-cont cataloged budget **plus** the uncataloged speckle companion), and what candidate radius does it imply under each defensible stellar-radius scenario? What did QLP actually correct for in the reported 2910 ppm?

### Multi-epoch consistency question — CORE / HIGH-VALUE BONUS

Are depth, duration, shape, and timing mutually consistent across Sectors 14, 41, 54, and 81 (2019–2024) after per-sector cadence, aperture, and dilution differences are modeled? The required answer uses predeclared per-sector estimates and an injection-derived scatter yardstick; if the method ladder promotes a hierarchical model, does its between-epoch variance posterior reach the same conclusion?

### GMU-night question — CORE (program requirement)

Does the correctly identified GMU R-band time series show any eclipse-like event in its phase 0.249–0.331 window; what upper limit does it place; what photometric precision was achieved (benchmarked against simultaneous Sector 54 photometry); and why was the night scheduled off-transit (ephemeris archaeology)?

### Source-separation question (the 0.52″ pair) — HIGH-VALUE BONUS / ADVANCED

Do the public multi-band depths (MuSCAT2 g,r,i,z_s; CMO g′; OAA Ic; KeplerCam i′; ULMT r′; GMU R; TESS) show the chromatic signature expected if the eclipse is on the ΔI = 1.7 companion rather than the primary? Does the transit-derived stellar density, fit per stellar scenario, break the tie?

### Stellar characterization question — HIGH-VALUE BONUS

Can public data (Gaia DR3 astrometry with its RUWE 3.69 caveat, SED photometry with the blend modeled, TRES parameters, TESS rotation) decide between the ~1.34 R☉ dwarf and ~2+ R☉ subgiant scenarios — and if not, can the paper at least *bound* the scenarios honestly?

### Variability and time-series question — HIGH-VALUE BONUS

Does the TESS photometry show the rotation period predicted by v sin i (≲ 5–8 d), γ Dor/δ Sct-type pulsation, or flares; do any of these bias per-sector depths; and are transit windows free of variability artifacts?

### Ephemeris question — HIGH-VALUE BONUS

Do individual transit times across 2019–2024 support a linear ephemeris; what is the refined prediction covariance; and what are the concrete 2026 re-observation windows (screening arithmetic already gives ~2026-07-15 00:52 UT, ~2026-08-06 08:35 UT, ~2026-08-15 02:28 UT, σ ≈ 7–9 min — verify before use)?

### False-positive question — ADVANCED

How do the localization results, dilution budget, chromatic constraints, density consistency, imaging contrast curves, TRES multi-epoch RV stability, and (if shared) TFOP notes change the relative probabilities of planet-on-primary, planet-on-companion, EB-on-companion, hierarchical EB, background EB, and systematic scenarios — via TRICERATOPS(+) with the real contrast curves, plus a transparent scenario matrix?

### Physical-context question — CONDITIONAL

At P = 2.915 d, the candidate sits interior to the published Neptunian ridge (3.2–5.7 d), i.e., in the hot-Neptune desert — *if* the Neptune-scale radius survives dilution and stellar-radius scrutiny. If instead Rp inflates to Jupiter scale, it is an ordinary hot-Jupiter candidate. Which is it, and what follow-up decides?

## Hypotheses and outcome-neutral tests

The hypotheses must be written before fitting the final models — in particular, before inspecting per-sector depth posteriors side by side and before unblinding the GMU window fit.

- **H1:** A separate fit to each TESS sector recovers a transit at the shared period with SNR sufficient for a per-sector depth estimate.
- **H2:** Per-sector dilution-corrected depths are consistent within the predeclared injection-derived scatter. If a hierarchical model is promoted, its between-sector variance parameter should support the same verdict.
- **H3:** Difference-image/centroid localization in each sector is consistent with the target pair's position, and no resolvable neighbor shows the signal in the SG1 NEB data or our reanalysis of it.
- **H4:** The GMU night (phases 0.249–0.331) is flat: no eclipse-like event above the frozen detection threshold; the achieved precision makes that limit meaningful.
- **H5:** The multi-band depth set is achromatic within uncertainties (planet-on-primary expectation), as opposed to the chromatic slope predicted under companion-hosted scenarios (compute the predicted slopes *before* measuring).
- **H6:** The transit-derived stellar density is consistent with at least one stellar scenario; a density consistent with *no* scenario indicates blend, eccentricity, or model failure.
- **H7:** A single timing offset is not evidence of TTVs; a TTV interpretation requires multiple precise timings and model comparison.

### Recommended primary endpoints and classification rule

Two frozen primary endpoints, reflecting the two data domains:

- **TESS endpoint:** the predeclared per-sector depth-consistency statistic compared with the injection-derived null yardstick, plus the per-sector localization verdicts. If a hierarchical model is promoted, its between-sector variance posterior is an additional endpoint, not a replacement for the transparent comparison. "Consistent" requires agreement under the applicable depth test and no sector's difference image excluding the target pair.
- **GMU endpoint:** the posterior upper limit (declared credibility level) on any eclipse depth within the observed window, from the frozen transit-plus-baseline model with duration fixed to the TESS-derived value, plus the achieved out-of-transit RMS benchmarked against simultaneous S54 photometry.

The mentor should approve exact thresholds before unblinding. Strong default framework:

- **Consistent multi-epoch signal:** all four usable sectors recover the signal; depth scatter is compatible with the injection yardstick (and hierarchical variance is compatible if that model is promoted); localization verdicts are uniform; density is consistent with ≥1 stellar scenario.
- **Tentative consistency:** one sector weak or one diagnostic (localization, chromaticity, density) marginal.
- **Physically inconsistent signal:** between-epoch variance excludes zero robustly, or localization/chromaticity/density excludes the clean interpretation — this is a *positive result* (false-positive evidence), not a failure.
- **Inconclusive:** data quality prevents the tests from discriminating.

Do not combine many weak, correlated diagnostics as though they were independent sigma evidence. Report each component and the classification logic.

## Claim ladder

Every statement in the paper should stop at the highest rung supported by the completed evidence.

| Rung | Permitted claim | Minimum evidence |
|---|---|---|
| 0 | No interpretable constraint | Wrong target, insufficient baseline, severe systematics, or inadequate precision. |
| 1 | A recurring dip exists in TESS data | Descriptive; already established by the TOI alert — reproducing it is not a contribution by itself. |
| 2 | The signal is separately recovered per epoch | Own extraction and fits in each of the four sectors with defensible uncertainties. |
| 3 | The signal is physically self-consistent across epochs | Hierarchical consistency posterior + duration/density coherence under a declared dilution model. |
| 4 | The signal is localized to the unresolved TIC 390988385 pair | Difference images/centroids across sectors + NEB re-derivations exclude resolvable neighbors to stated limits. **This rung, not rung 5, is the realistic ceiling of seeing-limited photometry for this target.** |
| 5 | The signal is attributed to the primary (or the companion) | Chromatic-depth analysis, density consistency, AO-resolved photometry, or RV/TFOP evidence actively discriminating between the pair members. |
| 6 | Statistically validated planet | Full validated-planet criteria including the companion scenario quantitatively retired; vetted inputs; mentor approval; no contradictory evidence. Likely unreachable with current public data — say so rather than stretch. |
| 7 | Dynamically confirmed planet | Mass or dynamical evidence meeting accepted standards; the unpublished HIRES RVs could do this, but they are not ours — outside this paper unless shared and permitted. |

Rung 5 is already partially populated by the public imaging (companion census) and TRES (no obvious SB2 reported — verify from the actual spectra/plots, not the absence of a note). The paper's job is to add genuinely distinct measurement layers or useful implementation stress tests, while labeling which is which.

## Analysis governance and anti-bias plan

### Freeze before unblinding

Before fitting the GMU window or comparing per-sector depth posteriors, save a dated configuration containing:

- Adopted coordinates and ephemeris source (and the reconstructed scheduling-era ephemeris, separately).
- Dilution model: cataloged contamination treatment, companion flux ratio and its band-dependence, and how each pipeline's prior corrections are undone/reapplied.
- Calibration rules for the GMU night (including the `-final` flats decision).
- Initial aperture and background-annulus grids (GMU and TESS pixel masks).
- Comparison-star eligibility rules.
- Quality metrics and rejection thresholds.
- Primary detrending covariates (including a meridian-crossing breakpoint term for the GMU night).
- Primary transit model, priors, and the hierarchical-model structure.
- Detection/upper-limit and consistency criteria.
- Injection–recovery criteria.
- Primary figures and tables.

If practical, tune systematics models on the GMU night's first and last hours (guaranteed out-of-transit under both candidate ephemerides) and on out-of-transit TESS segments. Any change after unblinding belongs in the analysis ledger with its reason and effect.

### Separate confirmatory and exploratory results

- **Confirmatory:** the frozen primary reduction and models.
- **Robustness:** predeclared alternative reasonable reductions.
- **Exploratory:** ideas added after seeing results.

Exploratory findings may motivate future work but must not be presented as preplanned independent evidence.

### Preserve all forks

Do not keep only the aperture, comparison ensemble, dilution treatment, or stellar scenario that produces the most planet-like answer. Save the full analysis grid, including null and contradictory outcomes. For this target the temptation will be strongest in the dilution/stellar choices — which is precisely where the paper's credibility lives.

### Claim–evidence ledger — signature research-engineering contribution

Create one machine-readable row per atomic claim before drafting prose. The minimum schema is:

| Field | Meaning |
|---|---|
| `claim_id` | Stable identifier used in text, figures, and review notes |
| `claim_text` | One checkable statement, not a paragraph |
| `status` | Proposed / supported / contradicted / silent / unusable / superseded |
| `claim_rung` | Highest rung the claim would require |
| `observation_id` | Raw observation or catalog snapshot that bears on it |
| `measurement_channel` | GMU ground photometry, TESS photometry, high-resolution imaging, spectroscopy, catalog inference, or literature |
| `implementation` | AIJ, custom Python, QLP, SPOC, custom FFI extraction, manual catalog query, etc. |
| `dependence_notes` | Shared photons, shared catalog fields, inherited corrections, simultaneous observation, or other non-independence |
| `config_and_commit` | Exact configuration plus code revision |
| `quantitative_result` | Estimate, uncertainty/limit, units, and comparison rule |
| `limitation` | What the evidence cannot establish |
| `review` | Student verification date and mentor review status |

The ledger prevents three common errors: turning multiple reductions of the same photons into multiple independent confirmations, silently inheriting a catalog value through several pipelines, and upgrading a qualitative consistency check into a validation claim. The scenario evidence matrix later in this plan should be generated from this ledger rather than maintained as a separate hand-written truth.

### Measurement-dependence map

Predeclare these relationships before combining evidence:

- AIJ and custom Python reductions of the 2022 GMU frames are **independent implementations of one observation**, not two observations.
- QLP, SPOC/TESS-SPOC, eleanor, TGLC, and custom apertures may provide useful pipeline-ablation evidence, but extractions of the same TESS pixels do not multiply the astrophysical evidence.
- Sector 54 and the GMU night are distinct instruments and measurement processes observing simultaneously, but their astrophysical signals are correlated because they view the same star at the same time.
- TIC values may inherit Gaia or other catalog inputs; quoting both is not automatically catalog-independent corroboration.
- The public SG1 nights, high-resolution images, and spectra can be distinct observations, but reductions or summaries derived from one posted file must stay linked to that parent observation.
- A pipeline's centroid, crowding correction, and transit depth can share upstream pixels and model assumptions. Report agreement as a cross-check, not as three independent votes.

Use the word **independent** only for the dimension that is actually independent: observation, instrument, implementation, or analyst. If no data product is comparable to a claim, mark it silent or unusable; absence of a measurement is not contradiction.

### Double-entry verification for load-bearing values

Verify each load-bearing identifier, time, ephemeris value, frame count, dilution term, and reported depth by two paths—for example, header parsing plus manual FITS inspection, an archive API export plus the rendered target page, or a scripted calculation plus a hand-worked test case. Record and reconcile disagreements before modeling. This is intentionally limited to load-bearing values; duplicating every cosmetic metadata field would add work without adding reliability.

### Ablation plan

The robustness grid should answer which evidence stream or decision actually carries the conclusion:

1. drop one comparison star at a time and then the full comparison-star subgroup;
2. drop each discretionary quality cut;
3. drop each detrending covariate;
4. remove catalog contamination and unresolved-companion dilution separately;
5. fit each sector alone, then leave one sector out of the combined result;
6. compare official and custom apertures while labeling them as shared-data pipeline ablations;
7. remove each public SG1 observation from the chromatic synthesis;
8. repeat the classification using only genuinely distinct measurement channels.

Report the direction and magnitude of the change, including null changes. Do not rerun only the ablations that make the preferred interpretation look stronger.

## Complete data inventory

### Required ground-based inputs — CORE

- TOI-3505.01 GMU science FITS sequence for 2022-07-21/22 (283 frames — delivered; verify against any observatory log).
- Same-night R flats (10, delivered as `-final` — establish provenance), 50 s darks (10), 3.5 s flat-darks (10), focus frames.
- Observer log, weather notes, and clock-sync information for the night (Omegalambda logs, if retained — ask).
- Telescope, detector, gain, read-noise, pixel-scale, binning, and filter metadata (headers delivered; read noise not in headers — get the STX-16803 value and verify on darks).
- Finder chart and the sign-up spreadsheet row (scheduling-era ephemeris evidence).
- The 2021-06-28 GMU night: ExoFOP posted products (public), and the raw frames from program drives if they exist.

### TESS inputs — CORE for the thesis

- All available light curves and target pixel files for **Sectors 14, 41, 54, 81**: QLP, TESS-SPOC, and any SPOC 2-min products (enumerate with a recorded `lightkurve` search; do not assume cadences).
- **DV reports, mini-reports, and DV time series for the three TCEs** (`s0014-s0086`, `s0054`, `s0081`), with the producing pipeline identified for each.
- TESScut FFI cutouts for custom extraction cross-checks in all four sectors.
- Data-release notes for sectors 14, 41, 54, 81 (verify sector date ranges and any anomalies — S14's early-mission systematics and each sector's scattered-light windows matter in the galactic plane).
- TIC entries and contamination-related fields (`CROWDSAP`, `FLFRCSAP`, TIC contamination ratio 0.547163) for every product used.
- Cadence-level quality flags, centroids, background, and momentum-dump timings per sector.

### Catalog and follow-up inputs — HIGH-VALUE BONUS

- Current NASA Exoplanet Archive TOI snapshot (dated) and `pscomppars` null-result record.
- ExoFOP metadata and the ~147 files: all six SG1 photometry sets (measurement tables, notes, NEB tars), ShARCS images/contrast data, SOAR sensitivity curve, TRES products.
- Gaia DR3 cone search exports (queries recorded; the 2.5′ and 30″ queries used for this document are re-runnable verbatim) including RUWE, `non_single_star`, variability flags.
- TIC neighbor table with TESS magnitudes for the NEB/dilution arithmetic.
- RAVEN released tables (Zenodo 19661443) — membership check for TIC 390988385.

### Optional external inputs — CONDITIONAL

- Pan-STARRS/DSS/2MASS archival images (PM here is only ~3–5 mas/yr, so the background-star displacement test is weak — note it rather than oversell it).
- ASAS-SN and ZTF time series (G ≈ 11.5 is near/inside ZTF saturation — check flags) for long-term variability.
- Public SED photometry + Gaia DR3 XP spectrum for stellar characterization **with the 0.52″ blend explicitly modeled** — a naive single-star SED fit of blended photometry is worse than none.
- Upcoming TESS coverage via WTV/tess-point once the tool responds.

External surveys must be used within their licenses and cited. Restricted TFOP products require permission.

## Ground-based reduction — CORE

### 1. Raw-file audit

1. Build the immutable manifest (checksums before/after zip extraction).
2. Inspect header consistency and exposure chronology (283 frames, ~73.7 s cadence — flag any gap > 2 cadences).
3. Plate solve representative early, middle, and late frames — and frames bracketing the meridian passage (~frames 140–210) to quantify any field rotation/flip.
4. Confirm the target remains in the field and quantify drift.
5. Plot median background, FWHM, ellipticity, airmass, centroid position, and total counts versus time; measure peak ADU of the target and bright comparisons on early frames (saturation check at 50 s).
6. Mark clouds, focus changes, meridian effects, tracking jumps, saturation, cosmic rays, and edge proximity without looking at any light-curve shape.

### 2. Calibration audit

1. Inspect every calibration frame; reject only with a recorded reason.
2. Determine what `-final` means for the delivered flats (already dark-subtracted? combined-then-split? ask); document the adopted algebra: darks at 50 s for science, darks at 3.5 s for flats, no separate bias — this is a coherent scheme only if applied exactly.
3. Create master dark and flat products with robust combination; normalize flats and inspect gradients, dust features, and stability.
4. Confirm dark exposure and detector temperature compatibility (headers: −20.4 °C setpoint).
5. Compare calibrated and uncalibrated statistics for a frame sample; test whether flat-fielding introduces structured artifacts.
6. Save master products and processing logs.

### 3. Time standard

1. Preserve original DATE-OBS (exposure start, per header comment) and the header JD/JD_UTC/BJD_TDB/HJD columns as provenance.
2. Re-derive mid-exposure BJD_TDB independently with astropy from DATE-OBS + 25 s + site coordinates + target coordinates; reconcile with the header BJD_TDB (spot checks in this research phase agree to the level expected of a mid-exposure convention, but verify across the whole night, and establish the acquisition computer's clock-sync provenance).
3. Fit all models in one declared standard (BJD_TDB); retain UTC/JD columns for audit.

A one-minute time error is fatal to the ephemeris-archaeology argument; treat timing as a first-class deliverable.

### 4. AstroImageJ primary reduction

AstroImageJ remains the program-aligned primary workflow unless the mentor directs otherwise. Record:

- AstroImageJ version and operating system; calibration options and master-frame paths; plate-solve method.
- Target and comparison-star coordinates; aperture and annulus radii; centroiding; variable-aperture settings.
- Comparison ensemble and weighting; removed frames with reasons; detrending columns.
- Exported time, flux, uncertainty, and quality columns.
- Export the `.plotcfg` and measurement table for every Discord post so others can reproduce the plot (program requirement).

Guides: [AstroImageJ documentation](https://astroimagej.com/guides/legacy/) and Conti's *A Practical Guide to Exoplanet Observing* (ask the mentor for the current TFOP SG1 guidelines version).

### 5. AIJ NEB analysis — CORE for localization (honestly scoped)

Run the AstroImageJ NEB analysis exactly as TFOP SG1 does — light curves for every Gaia/TIC star within ~2.5′ bright enough to matter, with per-star clearance arithmetic — **while stating plainly that an out-of-transit night cannot clear neighbors of the transit signal itself**. What our night's NEB run delivers: neighbor variability at phases 0.25–0.33, deep-eclipse detections at alias periods, and a methods rehearsal. The *transit-relevant* NEB evidence comes from re-deriving the ULMT (2021-10-15) and KeplerCam (2022-06-14) NEB checks from their posted measurement tables. Take the mentors up on NEB help at office hours; this remains squarely inside program scope.

### 6. Independent Python cross-check

A separate, transparent Python extraction (astropy + photutils, already in `requirements.txt`) should reproduce the AIJ result without copying its outputs: FITS ingestion, source matching, aperture photometry with a declared background estimator, ensemble differential photometry, independent BJD_TDB, quality table, raw and detrended light curves, and the frozen window fit. Agreement strengthens confidence in the reduction but is not a second astrophysical observation; disagreement is a result to diagnose, not to hide. As an optional third implementation, EXOTIC is student-oriented and citable, but add it only if AIJ/custom disagreement remains unresolved—not to inflate a pipeline count.

## Differential-photometry design

### Comparison-star eligibility

Predeclare bounds based on the actual field (the galactic plane helps here — comparisons are plentiful):

- Unsaturated in every retained frame; adequate SNR.
- Far from bad pixels, edges, blends (check each candidate against Gaia for sub-seeing companions — this field is full of them), and strong gradients.
- Stable centroid; no known variability (cross-check Gaia variability flags and ASAS-SN/ZTF where usable).
- Similar brightness and color to the target when possible (BP−RP = 0.745), reducing differential extinction over the airmass 1.07–1.40 run.
- Present throughout the full sequence.

### Ensemble construction

Evaluate: equal-weight stable ensemble; inverse-variance weighting; iteratively reweighted using declared out-of-transit segments; leave-one-out jackknife; every comparison as a pseudo-target. Select the primary ensemble on stability and field properties, never on the flatness of the target curve.

### Differential flux

`f_rel(t) = F_target(t) / F_ref(t)`, normalized on a declared baseline. Propagate target, sky, read, and ensemble uncertainty rather than relying only on fitted scatter. Remember the target aperture necessarily contains the 0.52″ companion and, depending on radius, may clip the 5.7–8.4″ neighbors — tabulate exactly which Gaia sources fall inside each tested aperture.

## Full robustness matrix — HIGH-VALUE BONUS

Run and save a grid that changes one defensible decision at a time, then a limited set of combined alternatives.

### Aperture tests (GMU and TESS)

- Fixed apertures spanning the useful range; FWHM-scaled apertures (~0.8–2.5×) after field testing.
- Alternative background annuli; median vs robust-mode background.
- Centroid-tracked vs fixed-coordinate apertures.
- **TESS: recovered depth vs pixel-mask size is the single most diagnostic robustness plot for a 55%-contaminated field** — a depth that changes systematically as pixels are added is a dilution-model failure or a blend signature. Produce it for all four sectors.

### Comparison-star tests

- Primary ensemble; leave-one-out; bright-only and color-matched subsets; equal vs inverse-variance weights; single-comparison diagnostics; pseudo-target curves for every comparison.

### Quality-cut tests

- No discretionary cuts beyond invalid/saturated frames; frozen thresholds; slightly stricter and looser thresholds; contiguity checks (does any conclusion depend on one cluster of frames, e.g., around the meridian passage?).

### Detrending tests

- Constant, linear, quadratic time baselines; airmass; centroid x/y; FWHM; sky; total comparison flux; a meridian breakpoint term; limited physically motivated combinations; masked vs simultaneous fits; GP alternative (declared kernel/priors) as ADVANCED.

Do not include regressors merely because they flatten (or reveal) a dip. Compare models on out-of-transit predictive performance, information criteria where appropriate, residual diagnostics, and injection preservation.

### Robustness deliverable

A specification plot: fitted eclipse-depth limit (GMU) and per-sector depth (TESS) across every acceptable configuration, with the primary marked. Report the primary result plus the range across reasonable alternatives.

## Noise and uncertainty model — HIGH-VALUE BONUS

### Theoretical noise budget

Estimate target photon noise, ensemble photon noise, sky, read noise, dark-current uncertainty, flat-field uncertainty, and scintillation (Young 1967 as modified by Osborn et al. 2015; D = 0.813 m, t = 50 s, altitude 154 m, airmass 1.07–1.40 — screening estimate is a few×10⁻⁴ per frame, i.e., scintillation-competitive with photon noise for an 11th-mag star; compute properly). Compare with observed out-of-transit scatter; a large excess identifies unresolved systematics.

### Empirical white and red noise

- Unbinned residual RMS; RMS vs bin size against the white-noise line; β factor on ingress/duration timescales; autocorrelation; Allan-style stability; residual permutation; block bootstrap with correlation-motivated block lengths.
- **Unique lever for this night: simultaneous Sector 54 photometry.** Bin the GMU curve to the S54 cadence, compare variability and trends directly, and use TESS as an external reference for what the star actually did that night. Ground–space disagreement = ground systematics (or resolved-vs-unresolved flux differences — think before concluding).

Inflate parameter uncertainties when red noise is present. Do not quote a white-noise MCMC interval as the full uncertainty if the residuals are correlated.

### Null and negative-control tests

- Fit the eclipse model at many control times across the night; block-permute residuals; fit pseudo-targets at the same times; fit positive and negative boxes; run the full search on comparison stars to estimate false-alarm behavior.

## Statistical and time-series method ladder — CORE where simple, ADVANCED only when earned

The methodological spine is **verified inference**, not maximum model complexity. Use the least complex method that answers each scientific question, and promote a method only after its inputs and simpler benchmark are trustworthy. Every method below must satisfy the Lecture 1 understanding rule—if it cannot be explained at a whiteboard in office hours, it does not go in the paper. Statistics fundamentals and error propagation are covered in [Lecture 8](../data_and_lectures/Lecture8_Schar2025_clean_redesign.pdf); everything beyond it must be documented and validated on synthetic data first.

### Method promotion ladder

| Tier | Methods | Promotion rule |
|---|---|---|
| 1 — required | Descriptive QC, explicit uncertainty propagation, weighted least squares, robust scatter, phase folding, simple box/transit fits, bootstrap or permutation controls | Default paper methods; complete before advanced sampling |
| 2 — high value | Per-sector physical transit fits, injection–recovery, linear ephemeris with covariance, Lomb–Scargle/ACF variability checks | Use when synthetic recovery and residual checks pass |
| 3 — advanced | Hierarchical depth model, nested-sampling model comparison, GP/ARMA noise, scenario-level source posterior | Use only if it changes a named decision, survives prior/seed sensitivity, and can be explained and reproduced |

If Tier 1 and Tier 3 disagree, the result is not automatically "more accurate" because Tier 3 is sophisticated. Diagnose the discrepancy and prefer the claim that survives both.

### Bayesian components

1. **Per-sector transit fits** (`batman` forward model; emcee or dynesty sampling): depth, mid-time, duration/impact parameterization, limb darkening per band, dilution as an explicit parameter with a prior from the contamination budget (not a hard-coded correction), per-dataset jitter, finite-exposure integration (essential at 30-min S14 cadence).
2. **Hierarchical multi-epoch depth-consistency model:** per-sector true depths drawn from a population `depth_s ~ Normal(μ, τ)`; the posterior on τ *is* the consistency result. Partial pooling beats both "fit everything jointly and hide tension" and "fit separately and eyeball." Report prior-sensitivity for τ. Extend to the public SG1 depths (per instrument/band) once dilution and limb darkening per band are modeled.
3. **GMU window upper limit:** nested-sampling or Savage–Dickey model comparison between flat-plus-systematics and eclipse-plus-systematics; report the depth posterior's upper credible bound and the Bayes factor with its sensitivity to priors, not a bare p-value.
4. **Bayesian ephemeris fit:** linear (and, only if warranted, quadratic) ephemeris over all timings with full epoch–period covariance; posterior predictive transit windows for 2026.
5. **Scenario-level source assignment (ADVANCED):** posterior probabilities over {primary, companion, neighbor, systematic} given centroid, chromaticity, and density likelihoods — only with validated likelihoods; otherwise present the scenario matrix qualitatively.
6. **Inference hygiene everywhere:** multiple chains, R̂, effective sample size, recorded seeds, synthetic-data recovery before real-data claims, posterior predictive checks per dataset, deterministic-optimizer cross-checks, and covariance reporting (depth–dilution–impact-parameter–limb-darkening degeneracies are severe in exactly this regime).

### Time-series components

1. **Period search:** BLS and TLS per sector and combined; explicit alias/harmonic audit (P/2, 2P, 1-day aliases) — with a 55% contaminated aperture, alias confusion with a neighbor EB is a live scenario.
2. **Rotation/variability:** Lomb–Scargle and ACF per sector (v sin i predicts P_rot/sin i ≈ 4.6 d for 1.34 R☉, ~7.6 d for 2.2 R☉ — a detected rotation period is a quantitative arbiter between stellar scenarios); window-function and alias checks; sector-to-sector coherence; γ Dor/δ Sct frequency scan (Teff 6200–6400 K is in-range); flare inventory and transit-overlap audit.
3. **Red-noise modeling:** β factor, block bootstrap, and a declared GP (e.g., SHO kernel) or low-order ARMA as *diagnostics* of correlated noise; any GP used in the primary fit must have constrained, physically motivated hyperpriors.
4. **Change-point discipline for the GMU night:** the meridian passage is a predeclared candidate breakpoint; model it rather than discovering it post hoc.

### Guardrails

- No stacking of black boxes: every sampler, kernel, and prior is declared and defensible.
- No hierarchical model rescues bad photometry — the Bayesian layer sits *on top of* the frozen reduction, never in place of it.
- Bayes factors are prior-sensitive; report sensitivity ranges or do not report Bayes factors.
- A hierarchical model is not required for the paper to be strong. If four reliable per-sector depths are not available by the promotion deadline, use transparent per-sector estimates and an injection-calibrated consistency table.
- Never list a method in the abstract, activity description, or portfolio until its code, diagnostic output, and interpretation exist in the repository.

## TESS reanalysis — CORE for the thesis

### Four sectors, four data qualities

| Sector | Epoch | Camera/CCD | FFI cadence | Role |
|---|---|---|---|---|
| 14 | Jul–Aug 2019 | 1/3 | 30 min | Discovery-era baseline; the data behind the 2021 alert; worst cadence — finite-integration handling mandatory. |
| 41 | Jul–Aug 2021 | 1/3 | 10 min | First post-alert epoch. |
| 54 | Jul–Aug 2022 | 2/4 | 10 min | **Contains the GMU night** and the current TOI epoch; camera/CCD change vs S14/41 = different systematics and different contamination realization. |
| 81 | Jul–Aug 2024 | 2/3 | 200 s | Best FFI cadence; most recent epoch; extends the baseline to 5 years. |

(2-min/20-s target data: verify on MAST per sector; the s0054 and s0081 TCEs suggest SPOC processing but do not by themselves prove 2-min cadence. Record the answer with the search query.)

For each sector, separately:

1. Download all official light-curve products (QLP, TESS-SPOC, SPOC if present) and target pixel files/cutouts; record versions.
2. Read the data-release notes; map scattered-light and momentum-dump windows onto the transit times.
3. Compare SAP vs PDCSAP/KSPSAP and cross-pipeline extractions — **explicitly reconciling each pipeline's dilution treatment (`CROWDSAP`, QLP crowding) before comparing depths**.
4. Plot the full sector before masking anything.
5. Recover the period per sector with BLS/TLS; audit aliases.
6. Fit individual transits and the per-sector phase fold; feed the hierarchical model.
7. Record missing/partial events and gaps.

Then read the **three DV reports end to end** (s0014–s0086 multi-sector, s0054, s0081) and reproduce their headline diagnostics with a custom implementation: odd/even, secondary search, centroid offset, difference image. Pipeline-vs-own agreement is useful robustness evidence from shared observations; disagreement is a finding.

NASA guides: [TESS data products](https://heasarc.gsfc.nasa.gov/docs/tess/data-products.html) and [Target Pixel File tutorial](https://heasarc.gsfc.nasa.gov/docs/tess/Target-Pixel-File-Tutorial.html).

### Extraction sensitivity — the dilution centerpiece

- Official apertures; smaller/larger custom masks; difference-image-informed masks; background alternatives.
- Depth vs aperture size per sector (the blend-warning plot).
- Dilution correction on/off/varied: propagate the contamination budget's uncertainty into the depth posterior instead of applying a point correction.
- Alternative FFI extractions (eleanor, TGLC, `unpopular`) as pipeline-dependence checks — cite versions.

### Standard transit-vetting diagnostics

Odd/even depths and timings; secondary-eclipse search at phase 0.5 and elsewhere (a ΔI=1.7 companion EB scenario predicts potential secondaries — search hard); V-shape assessment (expected anyway at b ≈ 0.8 — do not over-read it, but do not ignore it); individual-event consistency; model-shift/uniqueness; in/out difference images; centroid motion; aperture-contamination assessment; alias checks; transit-correlated background/systematics search. Use the DV products and (optionally) DAVE-style diagnostics as pipeline cross-checks. When they reuse the same TESS observations, their agreement is implementation robustness rather than a new independent measurement. A tool flag is not a classification.

### Stellar variability and flares

Mask transits; compute rotation/pulsation periodograms per sector; compare across sectors and with ASAS-SN/ZTF if usable; inventory flares and overlap with transits; test whether local trends bias per-sector depths (this couples directly into the hierarchical model's between-sector variance — a variance excess *caused by* unmodeled spots is a different conclusion than one caused by a blended EB). Treat any rotation-based age or spin-orbit speculation as exploratory context only.

## Gaia neighborhood and blend analysis — CORE for the thesis

The Gaia DR3 cone searches (2.5′ and 30″, run live 2026-07-14) establish: the target is **Gaia DR3 1824561377891719424** (G = 11.476, BP−RP = 0.745, **RUWE = 3.694**, parallax 1.2448 mas, `non_single_star = 0`, no variability flag); the 0.517″ speckle companion has **no separate Gaia entry**; the nearest cataloged neighbors are 19th–20th magnitude at 5.7″+; the brightest source within 30″ after the target is G = 17.52 at 8.42″; and the 2.5′ field contains >1000 sources (galactic plane).

### Preliminary neighbor screen (screening arithmetic — not a validation result)

Assumes full contaminant throughput and Gaia G as a band proxy: `flux ratio = 10^(−0.4 ΔG)`; `required eclipse ≈ observed depth / flux ratio` (using the cataloged-contamination-corrected depth where appropriate).

| Source | Sep | ΔG (or ΔI) | Flux ratio | Eclipse needed to mimic 2910 ppm | Verdict at screening level |
|---|---:|---:|---:|---:|---|
| **Speckle companion** | **0.52″** | **ΔI = 1.7** | **0.209** | **≈ 1.7%** | **Fully plausible mimic. Unresolvable by any seeing-limited data. The scenario to beat.** |
| Neighbor | 5.70″ | 8.57 | 3.7e-4 | ~790% | Impossible. |
| Neighbor | 6.33″ | 7.73 | 8.1e-4 | ~360% | Impossible. |
| Neighbor | 6.88″ | 6.73 | 2.1e-3 | ~140% | Impossible. |
| Neighbor | 8.42″ | 6.04 | 3.8e-3 | ~76% | Formally possible via near-total EB; extreme — test via SG1 NEB re-derivation and TESS difference images. |
| Wider field (30″–2.5′) | — | — | — | — | Must be enumerated via TIC/TESS-cont, not hand arithmetic: the 0.547 contamination ratio means substantial flux enters from beyond 30″ through the TESS PRF. |

The final analysis must replace this with TESS magnitudes, measured PRF/aperture throughput per sector, bandpass-aware contrasts, and uncertainties. TESS pixels are ~21″; GMU resolves most cataloged neighbors but never the 0.52″ pair. See [NASA's TESS crowding guide](https://heasarc.gsfc.nasa.gov/docs/tess/UnderstandingCrowdingv2.html).

### Required localization products

- GMU field image with target, Gaia/TIC neighbors, apertures, annuli, and N/E arrows.
- Our-night NEB table (variability-scoped, as discussed) + re-derived ULMT and KeplerCam NEB verdicts.
- Light curve for every neighbor bright enough to matter, per night.
- TESS pixel maps with apertures and Gaia overlays for all four sectors (`tpfplotter`).
- Per-sector difference images and centroid offsets (own + DV).
- Depth vs aperture size (all four sectors + GMU).
- Explicit unresolved-separation limits: GMU seeing/plate scale; TESS PRF/centroid precision; then what the ShARCS and SOAR contrast curves exclude interior to those limits — and what they do *not* exclude (the 0.52″ ΔI=1.7 companion is a detection, not an exclusion).

### Advanced contamination extensions

- **TESS-cont decomposition** (per-neighbor flux contributions per sector aperture; version cited) — this upgrades the 0.547 point value into a source-by-source budget, and is near-mandatory here.
- Pixel-level eclipse injection into specific neighbors to measure leakage into the aperture.
- Bandpass conversion of the companion's ΔI/ΔJ/ΔKs into ΔT_TESS, ΔR, Δg′… with stated assumptions (the two-band AO Δmags give a crude companion color — use it).
- Bayesian source assignment across {primary, companion, 8.4″ neighbor, systematics} (ADVANCED; only with validated likelihoods).
- Gaia epoch-astrometry caveats: RUWE 3.69 means the DR3 parallax/PM are unreliable; check for a DR4-era update before submission; do not build any quantitative claim on the DR3 parallax without an inflation model.

## Stellar characterization stress test — HIGH-VALUE BONUS (new section, load-bearing)

The candidate radius is quoted as 7.46 R⊕, but that number silently assumes (a) the TIC dwarf-ish radius, (b) the DR2-era distance, (c) zero flux from the uncataloged companion, and (d) the QLP dilution treatment. Each assumption is currently contestable. Deliverables:

1. **Scenario table (screening values verified 2026-07-14; recompute properly):**

| Stellar scenario | Basis | R* | Implied Rp (undiluted) | With +21% companion dilution (planet on primary) | Interpretation |
|---|---|---:|---:|---:|---|
| TIC dwarf | TIC v8.2, d = 373 pc | 1.34 R☉ | ~7.5 R⊕ | ~8.2 R⊕ | Hot-Neptune-desert object — rare if real. |
| TRES subgiant | log g 3.875 ± 0.109, M ≈ 1.2–1.4 M☉ | ~2.2 R☉ | ~12 R⊕ | ~13 R⊕ | Ordinary hot Jupiter. |
| Gaia-parallax giant-ish | π = 1.245 mas taken at face value + Teff | ~3 R☉ | ~17 R⊕ | ~18 R⊕ | Inflated-Jupiter/EB-suspicion regime. |

2. **Transit-derived stellar density as arbiter:** the reported 2.004 h duration at P = 2.915 d implies b ≈ 0.78 for the TIC radius and near-grazing (b ≳ 0.9) for the subgiant — or a smaller host (the companion?) at lower b. Fit ρ* from the transits per scenario and compare with spectroscopic log g. Predeclare the comparison; it is one of the few tools that can *separate* the pair members without new data.
3. **Blend-aware SED:** fit the combined photometry (Gaia/2MASS/etc.) as a two-star model with the AO/speckle Δmags as priors — a single-star SED here is methodologically indefensible.
4. **Honest distance treatment:** present both distances, the RUWE, and the consequences; refuse to pick a winner unless the evidence actually picks one.
5. **Rotation cross-check:** a measured TESS P_rot combined with v sin i = 14.74 km/s bounds R* sin i — an independent, nearly free constraint on the scenario table.

## Transit modeling — CORE / ADVANCED

### Primary model

Physically integrated transit model ([`batman`](https://arxiv.org/abs/1507.08285)) with: mid-time; radius ratio; scaled semimajor axis + impact parameter (or duration parameterization); limb darkening per band (TESS, R, and each SG1 band used); **dilution as a free parameter per dataset with contamination-budget priors**; baseline/systematics terms; per-dataset jitter; finite exposure integration (mandatory at 30-min S14 cadence). Catalog values enter as priors with cited uncertainty, not as fixed truths.

### GMU-only fits

1. Frozen upper-limit fit (duration/shape priors from TESS; mid-time scanned across the window).
2. Systematics-only null model.
3. Sensitivity fits (limb darkening, breakpoint on/off).

### TESS-only fits

- Separate per-sector fits (feeding the hierarchical model if promoted); shared-geometry joint fit as comparison; individual mid-times for timing; per-sector depths *before* any shared-depth enforcement, with pipeline dilution treatments reconciled explicitly.

### Joint and cross-instrument fits

- Joint TESS + public-SG1 + GMU fit with shared physical parameters, per-band limb darkening and dilution, per-dataset noise/baselines.
- The chromatic test: alternative fit with per-band depths free → depth-vs-wavelength slope vs the companion-blend prediction.
- Posterior predictive checks per dataset.

### Inference safeguards

Identifiability and prior-sensitivity audits; multiple chains, convergence diagnostics, ESS, seeds; deterministic-optimizer cross-check; synthetic-data recovery; report depth–dilution–b–limb-darkening covariances; never infer mass, density, composition, or atmosphere from photometry alone.

### Program-aligned cross-checks

The paper-structure lecture names **ExoFAST** as an expected tool; an ExoFASTv2 (or `juliet`/`allesfitter`) cross-fit of the primary `batman` fit is both robustness and program alignment. Cite what is actually used; understand every setting.

## Ephemeris and timing analysis — HIGH-VALUE BONUS

### Planning propagation (verified arithmetic, 2026-07-14)

For cycle `N`: `T_N = T0 + N·P`; `σ(T_N) = sqrt[σ(T0)² + N²σ(P)²]` under independence.

- GMU-night window: phases 0.2485–0.3310; nearest transits at cycle −4 (2459781.8738) and −3 (2459784.7889); back-propagated σ ≈ 3.0 min. **No transit in the window.**
- Summer 2026 screening windows (verify with the Swarthmore transit finder + refined ephemeris before requesting anything): ≈ 2026-07-15 00:52 UT (cycle 495), ≈ 2026-08-06 08:35 UT (cycle 503), ≈ 2026-08-15 02:28 UT (cycle 506); σ ≈ 7–9 min at 1σ. Transits of this target are observable from GMU *during this program cycle* — see the gated extension.

### Timing workflow

1. Fit every usable TESS transit time across Sectors 14, 41, 54, 81 with a consistent shape model.
2. Assign integer epochs carefully across the 2019–2024 baseline (~630 cycles; a mis-count is a silent catastrophe — verify by residual continuity).
3. Weighted linear ephemeris with epoch–period covariance; O−C; leave-one-out influence.
4. The archaeology deliverable: S14-only ephemeris propagated to 2022-07-22 vs the current ephemeris — the figure that explains the GMU night.
5. Compare pre-S81 vs post-S81 prediction uncertainty for 2026 ("future-observability rescue" framing).
6. Quadratic term only if the model comparison genuinely demands it.

Student-led ephemeris refinement is a real contribution when the timing is trustworthy; see the [ORBYTS ephemeris-refinement study](https://arxiv.org/abs/2005.01684) and [TESS ephemeris-recovery analysis](https://arxiv.org/abs/1906.02197).

### TTV guardrail

A single nonzero O−C point is a timing offset, not a TTV detection. A TTV claim requires repeated coherent deviations and a model comparison against no-TTV.

## Injection–recovery — HIGH-VALUE BONUS

Two distinct campaigns, because the two data domains answer different questions:

### GMU night (upper-limit calibration)

- Inject physically shaped eclipses (0.5–10 ppt; durations around 2.0 h; mid-times across the window) into minimally processed data before the reduction choices under test; recover with the frozen pipeline.
- Deliverable: completeness vs depth at the covered phases → converts "we saw nothing" into "we exclude eclipses deeper than X at these phases with Y% completeness."

### TESS sectors (consistency calibration)

- Inject transits of known, fixed depth into each sector's pixel/flux data; measure per-sector recovered-depth bias and scatter under the full extraction+dilution pipeline.
- Deliverable: the expected between-sector depth scatter under the null hypothesis "one constant astrophysical signal + our pipeline" — the yardstick against which the hierarchical τ posterior is judged. Without this, a nonzero τ is uninterpretable.

Freeze recovery rules before running the final grids. If GMU completeness at 3 ppt is poor, the honest statement is "the night does not constrain eclipses at this depth," not a fabricated limit.

## Formal false-positive validation — ADVANCED

Attempt only after the reductions, localization, and dilution work are correct. Statistical validation is conditional on inputs and scenario models; it is not proof — and for this target it may be structurally out of reach.

### TRICERATOPS / TRICERATOPS+

Use the current maintained version; cite the release. The 2025 extension ingests ground-based light curves in separate bandpasses — exactly what the public SG1 archive provides; see [Giacalone et al. 2021](https://arxiv.org/abs/2002.00691) and [TRICERATOPS+ (Gomez Barrientos et al. 2025)](https://arxiv.org/abs/2508.02782).

Inputs largely exist publicly: TESS pixel data and apertures (four sectors); Gaia/TIC neighbor properties; own transit fits; **ShARCS J/Ks and SOAR I contrast curves; the detected 0.52″ companion as an explicit resolved-companion scenario** (do not let the tool treat the field as companion-free); SG1 multi-band light curves; TRES constraints via mentor/TFOP notes. Record versions, seeds, priors, draw counts; repeat runs for numerical stability; report FPP and NFPP with per-scenario contributions; never quote only the most favorable run.

Expectation-setting: with a detected, unresolved, ~21%-flux companion, the planet-on-companion and EB-on-companion scenarios will likely dominate the FPP unless the chromatic/density evidence actively suppresses them. If FPP lands above validation thresholds, that is the finding — report it as such. The published thresholds (FPP < 0.015, NFPP < 10⁻³) are method-specific criteria, not permission to use the word "validated."

Stability audit: seed scatter; input perturbations within uncertainty; alternative priors; constraint inclusion/exclusion effects; double-counting checks (the companion must not enter as both a contrast-curve exclusion and a resolved star); preserved scenario samples.

### `vespa` historical cross-check

Only as a transparent secondary calculation with [Morton et al. 2016](https://arxiv.org/abs/1605.02825) cited — unmaintained; treat as historical. Investigate disagreements with TRICERATOPS rather than averaging them.

### Scenario-level modeling

Explicitly compare: planet on primary; planet on the 0.52″ companion; EB on the companion; hierarchical EB; background/foreground EB (the 8.4″ G=17.5 star and the unenumerated wider field); planet/EB on a resolvable neighbor leaking through the aperture; instrumental/systematic artifact. State which observation excludes which scenario at what confidence, and which scenarios remain standing at the end. That closing table *is* the paper's contribution.

### Validation stop conditions

Do not claim validation if: identity is unresolved; the companion scenario is untested by any resolving observation; any difference image points away from the target; odd/even or secondary evidence favors an EB; FPP/NFPP is prior- or seed-unstable; required inputs are missing from the adopted standard; or the mentor has not approved the wording. Given the companion, plan the paper so it *does not need* the word "validated" to succeed.

## Physical and population context — CONDITIONAL

Strictly conditional on the dilution/stellar outcome, because the candidate's identity swings between categories:

- **If Rp ≈ 7–9 R⊕ survives:** the candidate sits at P = 2.915 d, *interior* to the published Neptunian ridge (3.2–5.7 d) — i.e., in the hot-Neptune desert ([Castro-González et al. 2024](https://arxiv.org/abs/2409.10517)). Desert residents are rare and disproportionately valuable; the metal-rich host and existing HIRES data sharpen the follow-up case. This outcome would make the candidate genuinely important — which is precisely why the skeptical analysis must come first.
- **If Rp inflates to ~12+ R⊕:** an ordinary hot-Jupiter candidate around an evolved F star; the interesting result becomes the *correction* of the catalog radius and the demonstration of how crowding + stale stellar parameters manufacture false desert candidates.
- **Either way:** produce the dated, reproducible period–radius figure with published boundaries, the candidate marked with scenario-dependent radius extent (an error bar spanning 7.5–18 R⊕ is itself an honest, striking figure), and selection-effect caveats.

Use predicted mass/RV (43.4 M⊕ / 17.2 m/s under the Neptune scenario) only to motivate future work, noting the prediction collapses under the Jupiter scenario (K would be ~10× larger — which the HIRES data would have trivially measured; their 4-year silence is itself weak evidence worth one careful sentence after mentor consultation).

## Complete extension catalog

| Extension | Level | Scientific value | Dependency | Concrete output |
|---|---|---|---|---|
| Claim–evidence ledger + dependence map | CORE | Prevents correlated inputs and shared reductions from inflating the conclusion | Frozen questions and manifests | Machine-readable ledger, channel map, generated scenario table |
| Sector 81 analysis + three-DV reproduction | CORE/BONUS | Newest epoch; pipeline vetting comparison | Public MAST products | Per-sector fits + DV-vs-own diagnostics table |
| TESS-cont contamination decomposition (×4 sectors) | HIGH-VALUE BONUS | Converts 0.547 into a source-by-source budget | Public tool + Gaia | Contamination decomposition figure/table |
| Companion-aware dilution model | CORE for thesis | The uncataloged 0.52″ star enters every depth | ShARCS/SOAR Δmags | Band-dependent dilution priors + corrected depths |
| Public SG1 multi-band depth synthesis | HIGH-VALUE BONUS | The chromatic blend test, zero observing cost | ExoFOP login + observer courtesy | Depth-vs-wavelength figure across 7 datasets |
| Re-derived ULMT + KeplerCam NEB checks | HIGH-VALUE BONUS | Transit-relevant localization our night cannot give | Posted measurement tables | Per-star clearance re-derivation |
| Prior GMU night (2021-06-28) reanalysis | HIGH-VALUE BONUS | Two-epoch in-house comparison; on-transit epoch | ExoFOP files; raw frames via mentor | Second GMU light curve + depth |
| Simultaneous S54-vs-GMU benchmark | HIGH-VALUE BONUS | Rare same-night space/ground cross-check | S54 products | Precision/variability comparison figure |
| Ephemeris archaeology (S14-only propagation) | HIGH-VALUE BONUS | Explains the missed transit; program-useful | S14 refit + spreadsheet row | Stale-vs-current prediction figure |
| Hierarchical multi-epoch depth model | HIGH-VALUE BONUS / ADVANCED | The consistency thesis, done properly | Per-sector fits + injection yardstick | τ posterior + specification plot |
| Rotation/pulsation characterization | HIGH-VALUE BONUS | Tests stellar scenarios; guards depth fits | TESS light curves | P_rot / pulsation report + R* sin i bound |
| Blend-aware two-star SED fit | ADVANCED | Additional stellar scenario constraint | Public photometry + AO Δmags | SED posterior under two-star model |
| Transit-derived ρ* vs spectroscopic log g | HIGH-VALUE BONUS | Separates pair-member scenarios cheaply | Own fits + TRES values | Density-consistency figure |
| RAVEN membership check | CONDITIONAL | Current-literature ML comparison if present | Zenodo tables | Ledger entry ± scenario-posterior comparison |
| TRICERATOPS+ with real companion | ADVANCED | Quantitative scenario probabilities | Complete inputs | FPP/NFPP + stability audit |
| `vespa` cross-check | ADVANCED | Historical framework comparison | TRICERATOPS done | Evidence matrix row |
| Gaia DR4 astrometry re-check | CONDITIONAL | RUWE 3.69 case may resolve | DR4 availability | Updated distance/multiplicity entry |
| ASAS-SN/ZTF long-baseline variability | CONDITIONAL | Rotation/EB context | Saturation check | Periodogram + phase curve |
| Archival image PM test | CONDITIONAL (weak here) | Background-star check | PM only ~3–5 mas/yr — low power | One appendix panel, honestly caveated |
| Alternative FFI extractions (eleanor/TGLC/unpopular) | HIGH-VALUE BONUS | Pipeline-dependence test | Same raw products | Cross-pipeline depth comparison |
| ExoFASTv2 / juliet cross-fit | HIGH-VALUE BONUS | Program-named modeling cross-check | Public tools | Parameter comparison table |
| GP systematics alternative | ADVANCED | Correlated-noise robustness | Constrained kernel | GP-vs-parametric comparison |
| BLS vs TLS + alias audit | HIGH-VALUE BONUS | Period integrity in a crowded aperture | TESS light curves | Recovery/alias plot |
| Individual-transit depth/time variation | CONDITIONAL | Consistency at event level | Enough adequate events | Event-by-event parameter plot |
| TTV search | CONDITIONAL | Dynamics | Multiple precise timings | O−C model comparison |
| Secondary-eclipse/phase-curve limits | HIGH-VALUE BONUS (EB test) | Constrains companion-EB scenarios directly | TESS precision | Eclipse-depth limit at phase 0.5 |
| Flare–transit interaction audit | CONDITIONAL | Prevents activity bias | Detectable flares | Flare mask + sensitivity result |
| Desert/ridge population figure | CONDITIONAL | Context, scenario-dependent | Dated archive query | Period–radius figure with scenario-spanning error bar |
| Reproducible research release | CORE for maximum value | Auditable and reusable; proves the methods actually exist | Permissions/licenses + tests | Code, manifests, configs, derived tables, tagged release |
| Exoplanet Watch contribution | CONDITIONAL | Community reuse | Mentor + data-owner permission | AAVSO submission |
| Educational methods appendix | HIGH-VALUE BONUS | Student-led reproducibility | Complete analysis log | Stepwise appendix |
| 2026 transit re-observation from GMU | GATED / FUTURE | Replaces the missed transit; live windows exist (Jul 15 / Aug 6 / Aug 15 2026, screening) | Mentor/telescope approval + refined ephemeris | New on-transit R-band light curve |
| Multi-color ground photometry (new) | GATED / FUTURE | Sharpens the chromatic test | Approved observations | Depth-vs-band point(s) |
| TFOP-note / HIRES-context integration | GATED | Potentially decisive (mass or FP verdict exists somewhere) | Mentor permission after Discord post | Restricted-constraint discussion, properly attributed |
| Rossiter–McLaughlin / atmosphere / JWST | GATED / SEPARATE | Future-work framing only | Confirmation + approvals | Future-work paragraph |

Machine learning has no load-bearing role in this project (the pivot away from the ML-pipeline thesis was correct—that framing belonged to the old target). Project AERIS already demonstrates environmental ML/LLM evaluation; duplicating that surface form here would make both projects less distinctive. The transferable contribution is the verification architecture: provenance, frozen decisions, dependence-aware evidence, ablations, negative controls, and claim ceilings. The only defensible near-term ML use is optional frame-quality triage with interpretable features and human review, plus someone else's RAVEN output as a documented comparison point if membership is confirmed. Neither may determine the scientific classification on its own.

## Outcome decision tree

| Observed outcome | Best scientific interpretation | Best next analysis | Claim ceiling |
|---|---|---|---|
| Four sectors consistent; localization clean; GMU flat as predicted; chromatic test achromatic | Signal is a stable, pair-localized transit candidate; primary-hosted scenario strengthened | TRICERATOPS+, density test, 2026 re-observation request | Rung 4–5: consistent, pair-localized candidate with primary-favoring evidence |
| Sectors consistent but chromatic slope detected | Companion-hosted scenario favored | Quantify implied companion eclipse; EB modeling | Nearby/companion false-positive evidence — a publishable positive result |
| Between-sector depths inconsistent beyond injection yardstick | Systematics, variability, or blend dynamics | Aperture/dilution audit; variability coupling; alias check | Physically inconsistent candidate — also publishable |
| A sector fails to recover the signal at expected SNR | Depth/duration wrong, or extraction issue | Injection into that sector; DV comparison | Diagnosed discrepancy |
| GMU night shows a significant dip at phases 0.25–0.33 | Off-ephemeris event: alias EB, neighbor variable, or systematics | Localize within night; alias periods; neighbor curves | Evidence against clean interpretation |
| Ephemeris archaeology shows S14-propagation predicted the night correctly | Scheduling was sound; miss was ephemeris drift | Quantified drift case study | Methodological contribution |
| Archaeology cannot explain the scheduling | Unresolved provenance | Document; ask mentor; report honestly | Descriptive only |
| Time metadata unreliable | Timing analyses invalid | Recover clock provenance or drop timing claims | Photometric shape only |
| A HIRES/confirmation paper appears mid-project | Novelty reframing required | Pivot to "independent public-data stress test vs published result" | Comparison study — still viable if the analysis is honest |

## Result-dependent title options

The working title — *Stress-Testing TOI-3505.01: A Dilution-Aware, Multi-Epoch Analysis Using TESS, Gaia, and Public TFOP Constraints* — is deliberately outcome-neutral. The final title still follows the achieved result, and the word "Validation" may appear only if the analysis reaches that rung with mentor approval. The program's guidance is to title simply and after the outcome ([Paperstructure, title guidance](../data_and_lectures/Paperstructure.pdf#page=5)). Outcome-matched variants:

### If multi-epoch consistency and pair-localization hold

**Stress-Testing TOI-3505.01: A Dilution-Aware, Multi-Epoch Analysis Using TESS, Gaia, and Public TFOP Constraints**

### If the companion/chromatic evidence dominates

**Source Ambiguity in a Crowded Field: Dilution-Aware Limits on the Origin of the TOI-3505.01 Signal**

### If inconsistency is the finding

**Multi-Epoch Tension in the Transit Signal of TOI-3505.01**

### If the ephemeris story leads

**Ground-Based Follow-Up of TOI-3505.01 and the Cost of Stale Ephemerides for TESS Candidate Vetting**

### Program-format fallback

**Ground-Based Light Curve Follow-Up Observations of TESS Object of Interest TOI-3505.01**

## Figure plan

The previous 19-figure main-text list would maximize activity, not communication. The paper should have **six main figures, seven only if the stellar result is genuinely decisive**. Use multi-panel figures when the panels answer one question.

### Main-text figures, in narrative order

1. **Observation geometry and timing:** full Sector 54 context, the GMU window, nearest current-ephemeris transits, and the scheduling-era/Sector-14 prediction if reconstructed.
2. **GMU observation and reduction:** plate-solved finder image plus raw/detrended differential light curve and the few diagnostics needed to justify retained data.
3. **Four-epoch TESS recovery:** per-sector phase folds with cadence and uncertainty made visible; no decorative combined fold that hides sector disagreement.
4. **Dilution and aperture dependence:** TESS/Gaia field map, source-by-source contamination budget, and recovered depth versus aperture size.
5. **Localization:** own difference-image/centroid results beside the relevant DV and re-derived NEB constraints, with the unresolved 0.52″ pair visibly marked as unresolved.
6. **Claim-deciding synthesis:** per-sector depths and injection yardstick; add the multi-band SG1 comparison or scenario evidence matrix in the same figure according to which one actually changes the conclusion.
7. **Optional stellar discriminator:** rotation/density/stellar-scenario result only if it separates plausible hosts or candidate radii.

### Appendix figures

Full-sector light curves; master calibrations; plate-solution overlays (including meridian-bracketing frames); GMU sky/FWHM/centroid/airmass panels; S54-vs-GMU benchmark; every pseudo-target curve; RMS-vs-bin-size; full model-selection and ablation grids; odd/even, secondary, model-shift diagnostics; individual transits; injection–recovery completeness; O−C diagram; contrast curves; MCMC convergence/posteriors; prior sensitivity; posterior predictive checks; null-time event-search distribution; scenario probability breakdown; period–radius context.

Every figure must state data source, time system, flux normalization, binning, excluded points, and raw-vs-detrended status. Export as PNG or vector formats, never JPEG (program rule).

## Table plan

Keep the main paper to five tables:

1. Target/host properties and public follow-up assets, with provenance and snapshot date.
2. GMU observing, calibration, timing, and frozen photometry configuration.
3. Neighbor census and dilution budget, including the unresolved companion and per-star mimicry arithmetic.
4. Per-sector/per-instrument results with uncertainties, data dependencies, and injection yardsticks.
5. Claim–evidence/scenario matrix generated from the ledger, including surviving explanations and the final claim ceiling.

Place the full file manifest, technical-event log, software versions, robustness grid, injection completeness, individual transit times, stellar scenarios, and sampler diagnostics in machine-readable supplements or appendices. A reader should be able to audit the result without forcing every audit record into the narrative.

## Paper structure aligned with the lectures

The program specifies the exact skeleton — Title; Authors and Affiliations; Abstract; Introduction; Observations; Analysis; Results; Discussion; Conclusions and Future Work; References — with connective-tissue paragraphs opening each section. Follow it literally; the mapping:

### Abstract

Four parts (context, aims, methods, results), written last, including the off-transit-night finding and whatever the consistency analysis actually shows.

### Introduction (five-paragraph formula)

1. Transiting exoplanets, TESS, and the candidate-vetting problem — with crowding/dilution as the featured failure mode.
2. Ground-based follow-up and TFOP; cite the program-provided example papers.
3. The need for this paper: TOI-3505.01's history (QLP faint-star alert 2021; parameters updated 2024), its five-year PC limbo despite a rich unpublished TFOP record, the crowded field, the uncataloged companion, and the absence of any published synthesis.
4. "In this paper, we present…" — the four-sector stress test, the GMU night, and the specific questions.
5. Section roadmap.

### Observations (data accounting only — no analysis or results)

- 2.1 TOI-3505.01 and host properties (TIC, Gaia with RUWE caveat, TRES, provenance-tagged).
- 2.2 TESS sector data (14, 41, 54, 81; cadences, pipelines, DV products, release notes).
- 2.3 GMU 0.8 m observations: 283 × 50 s R-band frames, 2022-07-22 02:13–08:00 UT, calibrations, site/instrument — a log, per the lecture checklist. Credit observers/automation per mentor guidance.
- 2.4 Public follow-up assets used (six SG1 photometry sets, ShARCS, SOAR, TRES; HIRES existence noted).

### Analysis (reproducible methods)

Calibration and time verification; differential photometry; NEB methodology (ours + re-derivations); TESS extractions and dilution model; variability analysis; transit and consistency methods at the highest promoted tier; robustness suite; injection–recovery; validation tools if used.

### Results (facts only)

Achieved precision (including S54 benchmark); the flat-window result and upper limit; ephemeris archaeology outcome; per-sector fits and the promoted consistency statistic/model; localization results; chromatic synthesis; density comparison; timing.

### Discussion

- 5.1 Interpretation: where the signal originates at each angular scale; which false-positive scenarios survive; what the candidate most plausibly is under each stellar scenario; the lecture's question list answered explicitly.
- 5.2 Context: desert-vs-hot-Jupiter framing conditional on radius; how the results interact with the unpublished TFOP record; comparison to other follow-up work; what one more observation (RVs already exist somewhere; a 2026 transit; AO time-series) would decide.

### Conclusions and Future Work

One-sentence answer per research question, quantitative evidence, claim-ladder ceiling, and the single most informative next observation.

### Writing-style checklist (from the paper-structure lecture)

Succinct, dry, active voice, "we"; consistent tense; every acronym defined at first use; no unjustified claims; no qualitative qualifiers or "meaningless mouthfuls"; connective tissue at every section boundary; legible captioned numbered figures, PNG/vector only; references as taught; AI-use disclosure per Lecture 1.

## Poster and talk plan (due August 1)

Build the poster from the frozen Phase 0–2 results, not from whatever is newest on July 30. Follow the paper's section flow on the program template (40×30″). Reserve real estate for: the GMU light curve with the ephemeris-archaeology annotation (this is the poster's story — a missed transit *explained* beats a routine detection for memorability), one TESS phase-fold panel, the pixel/contamination figure, and the multi-band depth or scenario-matrix figure. If selected for a talk: 5 minutes, ≤1–2 slides per section, practiced; talk selection depends on analysis posted to Discord — post early and often.

The spoken explanation should work at two depths: a 20-second answer ("the observation missed the transit, so I verified why and tested the candidate across four TESS epochs") and a technical answer that can defend timing, dilution, and evidence dependence. Do not mention applications, AERIS, or a personal brand on the scientific poster.

## Application-facing artifact strategy — INTERNAL, downstream of the science

The college value comes from a compact set of finished, verifiable artifacts—not from adding more analyses to an unfinished paper. Build these only from frozen results:

1. **Journal-style paper and symposium poster/talk:** primary evidence of domain learning and scientific communication.
2. **Public research repository:** a concise README, one-command or clearly staged reproduction path, tests, frozen configs, data lineage, and a tagged release. Do not publish restricted or observer-owned data without permission.
3. **Short technical case study (600–900 words):** explain one hard decision—preferably why agreement among pipelines is not the same as independent evidence, or how dilution changed the claim ceiling. Include one result, one failed/negative test, and one limitation.
4. **Two-minute walkthrough or three-slide explainer:** question → evidence architecture → bounded conclusion. Show the actual plots and repository, not an animated product demo.
5. **Mentor-verifiable contribution log:** dates, specific analysis/code ownership, questions brought to office hours, feedback incorporated, final deliverables, and any coauthor contributions. This supports accurate recommendations and activity descriptions without asking a mentor to reconstruct the summer later.

### Honest positioning rules

- Use the official program name and institutional relationship. Do not shorten the role to "NASA researcher" or claim employment by NASA unless the program explicitly authorizes that wording.
- Say "paper in preparation" only while a real manuscript exists; "submitted" only after submission; "published" only after public publication.
- Claim Bayesian, hierarchical, time-series, or machine-learning methods only when the relevant code, outputs, diagnostics, and interpretation are complete.
- Lead with the research question and contribution, not the prestige of the data source.
- Quantify inputs and artifacts (283 frames, four TESS sectors, tested pipelines, reproducible release) rather than inventing impact metrics for a candidate classification.
- Preserve the negative result. Discovering that the provided night did not cover the current transit is evidence of careful research, not something to hide.

### Portfolio narrative after completion

The strongest one-sentence bridge to the broader body of work is:

> I build verification systems for noisy physical-world data: in AERIS I test environmental explanations against process-distinct measurements, and in this project I test an exoplanet candidate across ground, satellite, imaging, and catalog evidence while tracking which sources are genuinely distinct.

Use this only as portfolio/application framing, not as paper language. The project should still stand on its own if AERIS is never mentioned.

## Proposed repository structure

```text
README.md                            # question, current status, result, quickstart, limitations
CITATION.cff                         # add for the tagged public release
data/
  raw/
    gmu/toi3505/20220721/          # extracted from the six delivered zips; manifest + checksums
    gmu/toi3505/20210628/          # if raw frames are recovered from program drives
    tess/toi3505/sector14/
    tess/toi3505/sector41/
    tess/toi3505/sector54/
    tess/toi3505/sector81/
  external/
    gaia/
    exoplanet_archive/
    exofop/                        # SG1 tables, NEB tars, ShARCS + SOAR products, TRES files
    dv_reports/                    # s0014-s0086, s0054, s0081
    raven/                         # if membership confirmed
  processed/
    gmu/toi3505/
    tess/toi3505/
  manifests/
    input_files.csv                 # source, checksum, license/permission, retrieval date
    catalog_snapshots/
    measurement_dependence.yaml
configs/
  toi3505_primary.yaml
  toi3505_robustness.yaml
  toi3505_dilution_budget.yaml
docs/
  toi-3505-research-maximization-plan.md
  claim-evidence-ledger.csv
  decision-log.md
  ai-assistance-ledger.md
  novelty-ledger.md
  public-case-study.md              # written only after results freeze
notebooks/
  01_data_identity.ipynb
  02_gmu_quality_control.ipynb
  03_ephemeris_archaeology.ipynb
  04_tess_sector_analysis.ipynb
  05_dilution_and_neighbors.ipynb
  06_transit_models_and_hierarchy.ipynb
  07_variability_timeseries.ipynb
  08_injection_recovery.ipynb
  09_false_positive_validation.ipynb
src/
  toi3505/
tests/
  unit/
  integration/
outputs/
  toi3505/
    figures/
    tables/
    models/
    logs/
    paper/
    poster/
    release/
```

Raw data remain immutable. The six delivered zips (~6.3 GB) currently sit in `data_and_lectures/` with Git LFS tracking configured in this repo — **decide deliberately** whether 6 GB belongs in LFS (quota, clone cost) or outside the repo with manifests and retrieval instructions; do not let the default win by inertia. Do not silently overwrite downloaded files.

The current repository has useful exploratory code and outputs but no top-level README, package structure, test suite, frozen configuration, or data-lineage manifest yet. Building that thin research-software layer is higher value than adding a fifth modeling framework. Keep notebooks for exploration; move claim-producing calculations and final figure generation into importable modules and scripts so they can be tested.

## Reproducibility standard

- Replace broad minimum-version ranges with a generated lock or exact environment export for the release; keep a readable top-level dependency list and record Python, OS, AstroImageJ, and Java versions.
- Record AstroImageJ version and every manual setting; archive the `.plotcfg` files posted to Discord.
- Save archive query text, retrieval date, and raw returned tables (the TAP/ExoFOP/Gaia/TESScut/exo.MAST queries used for this document are recorded here and re-runnable verbatim).
- Save random seeds and sampler settings; store machine-readable priors and configuration.
- Hash raw inputs, frozen configs, code commit, and claim-producing derived tables. A result freeze is valid only when those four identities are linked in a manifest.
- Paper figures generated from scripts, never hand-edited; plotting separate from numerics.
- Units and time standards in column names.
- Tests for: time conversion (header-BJD vs astropy-BJD), aperture flux on synthetic images, injection recovery, ephemeris cycle counting, and the dilution algebra (a sign error in a dilution correction is this project's most likely silent catastrophe — test it against hand-computed cases).
- Add one end-to-end integration test on a small, publishable fixture—not the 6 GB archive—that recreates a key table/figure summary.
- Run tests automatically on pushes if setup time permits; at minimum, record the exact test command and passing output at poster freeze and paper freeze.
- Clean-environment reproduction before submission; preserve negative/failed runs in logs; cite software and data releases; obtain permission before publishing GMU or restricted follow-up data.
- Tag the exact submitted state and create a release manifest. Archive only what licenses and program rules allow; a manifest plus retrieval instructions is preferable to republishing data without authority.

### Minimum public-release acceptance test

A new reader should be able to answer, from the repository alone: What was observed? Which files and catalog snapshots were used? Which decisions were frozen? Which outputs support the headline claim? Which evidence streams share upstream data? What failed or remained inconclusive? What command regenerates the main quantitative result? If any answer depends on private memory or a Discord thread, the release is not finished.

## Prioritized execution roadmap (calendar-mapped)

As of **July 20**, 12 days remain to the August 1 symposium and 40 days to the August 29 paper deadline. The old roadmap attempted nearly every plausible extension before either deadline. This revision protects a complete core paper, one distinctive verification contribution, and a small number of result-dependent bonuses.

### Phase 0 — Identity, provenance, and first public checkpoint (July 20–22)

- [x] Audit all six archives without silently treating duplicate filenames as new observations: 283 unique science frames, one duplicate science file, 10 × 50 s darks, 10 × 3.5 s darks, 10 flats, 11 unique focus FITS, and 22 focus PNGs.
- [x] Export the 283-frame header table and simple timing/airmass figures; verify that the current-ephemeris transit falls outside the observed window.
- [x] Prepare a plain-language Discord progress update and four simple images; **posting and mentor response remain pending**.
- [ ] Hash the six immutable archives and create `input_files.csv`; do not extract 6 GB into a second uncontrolled copy merely to satisfy a checklist.
- [ ] Plate solve early/middle/late and meridian-bracketing frames; confirm target identity against a Gaia/TIC overlay.
- [ ] Verify the calibration set and resolve the `-final` flats provenance.
- [ ] Create the frozen primary config, claim–evidence ledger, decision log, and measurement-dependence map.
- [ ] Post the current timing question/progress update in the program's light-curve or data-reduction channel, then record the answer instead of relying on memory.
- [ ] Export the 2026-07-14 Archive/ExoFOP/Gaia snapshot and queries into `manifests/`.

**Stop rule:** no interpretation until identity passes.

### Phase 1 — Internship light-curve core (July 21–25; use office hours)

- [ ] Calibrate the night in AstroImageJ; select apertures/comparisons on field quality; saturation check.
- [ ] Produce raw and detrended differential light curves; verify BJD_TDB independently.
- [ ] Compute and plot the night's phase coverage under current and reconstructed ephemerides; bring the archaeology question to the next available office hours.
- [ ] Run the NEB extraction on our night; download and re-derive the ULMT/KeplerCam NEB checks.
- [ ] **Post the actual light curve + settings + `.plotcfg` to Discord**—this, not the header-only checkpoint, unlocks TFOP context and supports talk eligibility.
- [ ] Reproduce the AIJ result in Python from the same frames; label this an implementation cross-check, not an independent observation.
- [ ] Add unit tests for ephemeris cycle counting, BJD convention, and one hand-computed dilution case.
- [ ] Draft outcome-neutral Observations and Analysis sections.

**Promotion rule:** if the light curve is not trustworthy by July 25, stop expanding the model stack and finish the ground-based methods/results cleanly with mentor help.

### Phase 2 — Poster-grade result (July 23–28)

- [ ] Download and manifest TESS products for Sectors 14, 41, 54, and 81 plus the three DV reports; settle cadence/product availability.
- [ ] Recover the signal separately in each usable sector with the simplest defensible fit; do not wait for a hierarchical model.
- [ ] Build the companion-aware dilution budget and one depth-vs-aperture plot, starting with Sector 54.
- [ ] Compare the GMU night with simultaneous Sector 54 photometry.
- [ ] GMU-window upper-limit fit + injection–recovery (GMU campaign).
- [ ] Fill the claim–evidence ledger for every poster statement and run the first leave-one-decision-out ablations.
- [ ] Choose the poster's claim ceiling from completed evidence; use "preliminary" where appropriate.

### Symposium checkpoint — poster/talk due August 1

- [ ] Freeze poster-input code, config, data manifest, and results by July 28.
- [ ] Build on the program template; practice the 20-second and 5-minute explanations.
- [ ] Run the minimum public-release acceptance test on the poster's headline result.

### Phase 3 — Paper core (August 2–12)

- [ ] Complete per-sector extraction, aperture/dilution tests, standard vetting diagnostics, and the relevant DV reproductions.
- [ ] Complete the GMU timing/precision result and ephemeris archaeology with its source provenance.
- [ ] Build the source-by-source contamination table and explicit unresolved-companion scenario.
- [ ] Run TESS injection tests sufficient to calibrate cross-sector depth comparisons.
- [ ] Generate the six main figures and five main tables from scripts.
- [ ] Draft Results as atomic statements linked to ledger rows.

### Phase 4 — One decisive bonus, not five (August 13–20)

Choose the first feasible option that changes the claim ceiling, then stop:

1. public SG1 multi-band depth/NEB synthesis;
2. hierarchical depth consistency, if promoted by the method ladder;
3. transit-density/rotation stellar discriminator;
4. TRICERATOPS+ only with complete, validated companion inputs.

Do not run `vespa`, a second sampler, a GP, and a two-star SED merely to lengthen a methods list. New observations remain gated future work unless the mentor identifies a specific feasible night that resolves the leading uncertainty.

### Phase 5 — Claim freeze, mentor audit, and submission (August 21–29)

- [ ] Freeze primary analysis and run the full tests/reproduction path.
- [ ] Result-matched title (drop "Validation" unless earned) and abstract; limitations and claim-ladder audit.
- [ ] Re-run novelty/catalog checks and obtain mentor review of the claim ceiling, role wording, authorship, and data-release permissions before Aug 25 if possible.
- [ ] Clean-environment reproduction.
- [ ] Send PDF to execed@gmu.edu and nasa.schar.program@gmail.com by Aug 29.
- [ ] Tag the submitted code state and preserve the exact submission PDF.

### Phase 6 — Public/application translation (after the scientific freeze)

- [ ] Write the short case study from the final ledger and one real methodological lesson.
- [ ] Finish the public README and permitted release; add `CITATION.cff`.
- [ ] Record the two-minute walkthrough or build the three-slide explainer.
- [ ] Give the mentor a concise contribution log and retain the approved role/project wording.
- [ ] Align the GitHub profile label with the actual project: astronomy / scientific computing / physical-data inference, not "geospatial."

## Value-versus-effort priority

### Must do

1. Correct data identity and calibration.
2. Trustworthy differential photometry and independently verified BJD_TDB.
3. The phase-coverage/ephemeris-archaeology result, stated plainly.
4. Per-sector TESS recovery with one declared dilution treatment and dependence-aware interpretation.
5. Claim–evidence ledger, frozen config, data manifest, and the three load-bearing unit tests.
6. Honest uncertainty, negative controls, and result-matched claims.
7. Discord posts early enough to obtain feedback and unlock program context.
8. A finished poster/talk and paper.

### Best scientific return for added work

1. Companion-aware dilution budget + TESS-cont decomposition.
2. NEB extraction (ours) + re-derived public NEB checks.
3. S54-simultaneity benchmark and depth-vs-aperture plots.
4. Multi-band depth synthesis from public SG1 data.
5. Injection-calibrated multi-epoch consistency; hierarchical only if promoted.
6. Ephemeris refinement 2019–2024 + 2026 windows.
7. Transit-derived density or rotation as a stellar-scenario discriminator.
8. Relevant DV diagnostic reproduction.

### Advanced only after the above is solid

1. TRICERATOPS+ scenario probabilities with complete companion inputs.
2. Hierarchical or Bayesian source assignment.
3. Blend-aware two-star SED fit.
4. GP alternatives, `vespa`, or TTV search.

### High-value downstream artifacts

1. Public reproducibility release and technical case study.
2. Mentor-verifiable contribution log.
3. Two-minute/three-slide explanation for a non-astronomy technical audience.

### Moonshots / separate projects

1. New on-transit GMU observation.
2. New AO time-series or multi-color campaign.
3. JWST track (separate paper, program-gated).

## Risk register

| Risk | Consequence | Early diagnostic | Mitigation |
|---|---|---|---|
| Wrong-field or mislabeled files | Study invalid | Header/plate-solve mismatch | Hard identity gate; Gaia overlay in a crowded field |
| Header BJD/clock provenance wrong | Archaeology argument collapses | astropy-vs-header reconciliation | Independent time derivation; treat headers as claims |
| `-final` flats hide unknown processing | Subtle photometric bias | Flat statistics vs raw expectations | Provenance question to mentor; test both calibration paths |
| Saturation of target/comparisons at 50 s | Nonlinear photometry | Peak-ADU audit on early frames | Aperture strategy; exclude saturated comparisons |
| Meridian-passage systematics | Fake dips/steps mid-night | Frame ~140–210 diagnostics | Predeclared breakpoint term; segment checks |
| Dilution algebra error or double-counting | Wrong depths everywhere | Hand-computed test cases | Unit-tested dilution module; one declared budget |
| Companion scenario untestable → validation impossible | Overclaim temptation | Claim-ladder audit | Plan the paper to succeed at rung 4–5; mentor wording review |
| Distance/stellar radius unresolved | Rp spans 7.5–18 R⊕ | RUWE, TRES-vs-TIC tension | Scenario table instead of false precision |
| Stellar variability (P_rot ~5 d) biases depths | Between-sector inconsistency misread | Rotation analysis before consistency verdict | Variability-coupled robustness tests |
| A HIRES/confirmation paper appears mid-project | Novelty framing wrong | Weekly ADS/arXiv alerts | Pivot to independent-comparison framing (decision tree) |
| Alias/neighbor-EB period confusion | Wrong physical model | BLS/TLS alias audit | Explicit harmonic checks; neighbor curves |
| Multiple testing across the robustness grid | Inflated significance | Analysis ledger | Frozen primary tests; null controls |
| Advanced model instability | Misleading FPP/posteriors | Seed/prior sensitivity | Repeat, cross-tool, simplify, disclose |
| Shared photons/catalog lineage counted as independent evidence | Confidence inflated without new information | Measurement-dependence map | Label observation vs implementation independence; channel-level ablation |
| Scope overload before Aug 1 / Aug 29 | Poster or paper unfinished | Roadmap status vs calendar | Finish phases in order; poster freeze July 28 |
| Application-driven scope creep | Scientific question gets replaced by a résumé narrative | Tasks added without changing a claim | Enforce the four-part task decision rule; translate only after freeze |
| AERIS language or an LLM layer forced into astronomy | Project reads as branding rather than domain research | Methods cannot be justified astronomically | Transfer verification habits only; keep AERIS out of paper/poster |
| Public claim outruns completed work or official affiliation | Credibility/ethics problem in applications | Activity wording names unrun methods or implies NASA employment | Evidence-backed wording audit with mentor; dated contribution log |
| Repository remains an exploratory dump | Work is difficult to verify or reuse | No README, tests, configs, or release path | Build the thin research-software layer before extra frameworks |
| 6 GB zips mishandled in git/LFS | Repo bloat or quota failure | LFS status check | Deliberate storage decision + manifests |
| Restricted TFOP data used improperly | Publication/ethics problem | Provenance audit | Permission or omission |
| TOI-3505/3506 conflation in program materials | Wrong-target contamination of searches | Double-check every identifier | TIC-first search discipline |

## Pre-mortem: the hostile referee test

Before submission, the draft must survive these attacks, each of which is currently live:

1. *"You analyzed a night with no transit and called it follow-up."* — Answer with the archaeology figure, the upper limit, the S54 benchmark, and the program's own outcome taxonomy. If those aren't in the paper, the referee is right.
2. *"Your localization claim ignores a 21%-flux companion at 0.5″."* — Answer with the claim ladder capped at the pair, plus whatever chromatic/density evidence actually discriminates. Never claim more.
3. *"Your depths are pipeline artifacts of inconsistent dilution corrections."* — Answer with the single declared budget applied uniformly, the depth-vs-aperture plots, and the injection yardstick.
4. *"Your stellar parameters are internally contradictory."* — Answer with the scenario table and the density test; do not present one radius as truth.
5. *"Your consistency claim is just wide error bars."* — Answer with the predeclared consistency statistic against the injection-calibrated null yardstick; include the τ posterior only if the hierarchical method was promoted.
6. *"This is all already known to TFOP."* — Answer with the novelty ledger: nothing is published, our night is undocumented, and no public synthesis exists. And verify that stays true the week of submission.
7. *"You counted three reductions of the same photons as three confirmations."* — Answer with the measurement-dependence map and channel-level ablations. If the prose still calls them independent measurements, the referee is right.
8. *"This is a software portfolio wrapped around a weak astronomy result."* — Answer with a paper whose question, figures, and conclusion are entirely astronomical; keep the software architecture in methods/reproducibility and the personal throughline outside the paper.
9. *"The public résumé says more than the repository proves."* — Answer with the tagged release, contribution log, and exact status language. Remove any method or outcome that cannot be opened and checked.

## Approval gates

| Gate | Approval/evidence needed | Work unlocked |
|---|---|---|
| Ground data identity | Conclusive identity audit | GMU reduction and interpretation |
| First light curve posted to Discord | Program/mentor review | TFOP context (incl. HIRES story); talk eligibility; full paper progression |
| Internal TFOP information | Explicit mentor permission | Restricted constraints or notes in the paper |
| ExoFOP file usage | Download + observer-courtesy/acknowledgment rules confirmed | SG1 tables, NEB tars, contrast curves in figures/validation |
| Statistical validation claim | Complete inputs incl. companion scenario + mentor review | "Validated" language, only if thresholds are genuinely met |
| New telescope observations | Mentor/telescope approval | 2026 on-transit night |
| JWST/atmosphere work | Lecture-defined approval after campus light curve | Separate project/paper |
| Public release | Data-owner and license permission | Repository data/products release |
| Public/activity wording | Completed artifact + mentor-verifiable role and status | Application, portfolio, and recommendation language |

## Definition of a successful project

The project succeeds if it produces a reproducible and correctly bounded answer to the localization and consistency questions — including the honest accounting of the off-transit GMU night. Success does not require a transit in our data, a validated planet, or a resolved companion.

A scientifically complete result includes:

- Proven data identity and auditable calibration/timing.
- The ephemeris-archaeology result, quantified.
- A GMU photometric benchmark with an injection-calibrated upper limit, cross-checked against simultaneous TESS data.
- Four separately executed sector analyses under one declared dilution model, with a consistency verdict calibrated by injections and a hierarchical result only if earned.
- Localization verdicts with explicit angular scales, including what is *not* excluded.
- A stellar scenario table replacing false precision; apply the density test if the data support it.
- A scenario evidence matrix stating what the candidate can and cannot be.
- Claims that stop at the supported rung and a named, justified next observation.

A maximum-value research-engineering result additionally includes:

- A frozen claim–evidence ledger and measurement-dependence map.
- Unit tests for the load-bearing timing and dilution arithmetic plus one small end-to-end reproduction fixture.
- Script-generated figures/tables, exact configs, input hashes, a tagged code state, and a clear public/private data boundary.
- A README and permitted public release that a new reader can audit without private memory.

The downstream application value is successful only if it can be expressed through completed artifacts and mentor-verifiable ownership. It should show one coherent strength—building trustworthy computational evidence for physical systems—while leaving the astronomy and environmental projects scientifically distinct.

## Immediate next actions

1. Hash and manifest the six existing zips; plate-solve representative frames and confirm identity.
2. Post the prepared header/timing checkpoint, then record the mentor response.
3. Resolve `-final` flat provenance and complete the first AstroImageJ reduction.
4. Create `claim-evidence-ledger.csv`, `measurement_dependence.yaml`, the decision log, and the frozen config before choosing a preferred light curve.
5. Reproduce BJD cycle counting and dilution arithmetic in tests.
6. Post the actual AIJ light curve with settings, measurement table, and `.plotcfg`; this is the program gate.
7. Reproduce that light curve in Python and diagnose any disagreement.
8. Download/manifest the four TESS sectors and three DV reports; start with the Sector 54 simultaneity and aperture test.
9. Audit the sign-up spreadsheet/program drives for the scheduling ephemeris and 2021 GMU night; request only the permissions/data actually needed.
10. Freeze the poster claim and inputs on July 28; postpone non-decisive advanced models.

## Reference starting set

### Program materials

- [Lecture 1: introduction and AI guidance](../data_and_lectures/Lecture1_intro.pdf)
- [Lecture 2: NASA missions](../data_and_lectures/lecture2_NASAmissions.pdf)
- [Lecture 3: working with data](../data_and_lectures/lecture3_workingwithdatapt1%20%281%29.pdf)
- [Lecture 4: light, telescopes, filters](../data_and_lectures/Lecture4_Schar2025_clean_redesign.pdf)
- [Lectures 5–6: telescopes, TESS, transits, false positives](../data_and_lectures/Lecture56_clean_redesign%20%281%29.pdf)
- [Lecture 8: research-project sequence, statistics, error propagation](../data_and_lectures/Lecture8_Schar2025_clean_redesign.pdf)
- [Paper structure](../data_and_lectures/Paperstructure.pdf)
- [GMU Schar Astro Scholars example papers](https://science.gmu.edu/academics/departments-units/physics-and-astronomy-department/observatory/schar-astro-scholars), including the [TOI-5372.01 example](https://science.gmu.edu/sites/default/files/2024-10/Ground-Based%20Follow-Up%20Observations%20of%20TESS%20Object%20of%20Interest%20%28TOI%29%205372.01.pdf) — and noting the JASR archive's TOI-3506.01 paper as a not-this-target caution
- Program-provided professional examples of ground-based TESS follow-up (Lecture 8, p.12): [high-school student example](https://emerginginvestigators.org/articles/ground-based-follow-up-observations-of-tess-exoplanet-candidates/pdf), [arXiv:2110.14344](https://arxiv.org/abs/2110.14344), [MNRAS 516, 4432](https://academic.oup.com/mnras/article/516/3/4432/6692879), [arXiv:2210.08179](https://arxiv.org/abs/2210.08179), [arXiv:2205.05709](https://arxiv.org/abs/2205.05709), [arXiv:2212.08242](https://arxiv.org/abs/2212.08242), [arXiv:2208.07328](https://arxiv.org/abs/2208.07328), and the mentor's TOI follow-up bibliography via [ADS](https://ui.adsabs.harvard.edu/search/q=%20author%3A%22Plavchan%22%20title%3A%22TOI%22&sort=date%20desc%2C%20bibcode%20desc&p_=0)

### Mission, catalog, and follow-up documentation

- [NASA Exoplanet Archive TOI table](https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=TOI)
- [ExoFOP-TESS target page for TIC 390988385](https://exofop.ipac.caltech.edu/tess/target.php?id=390988385)
- [TESS data products](https://heasarc.gsfc.nasa.gov/docs/tess/data-products.html) and [Target Pixel File tutorial](https://heasarc.gsfc.nasa.gov/docs/tess/Target-Pixel-File-Tutorial.html)
- TESS data-release notes for Sectors 14, 41, 54, 81 (retrieve and cite the specific DRN numbers)
- [Understanding TESS crowding](https://heasarc.gsfc.nasa.gov/docs/tess/UnderstandingCrowdingv2.html) and [TESS telescope information](https://heasarc.gsfc.nasa.gov/docs/tess/telescope_information.html)
- [TFOP ground-based follow-up overview](https://asd.gsfc.nasa.gov/archive/tess/ground_based_followup.html)
- [Gaia DR3 archive](https://gea.esac.esa.int/archive/), [documentation hub](https://www.cosmos.esa.int/web/gaia/dr3), and the Gaia DR3 RUWE/astrometric-quality documentation (needed to defend the RUWE 3.69 argument)
- [AstroImageJ documentation](https://astroimagej.com/guides/legacy/)
- [Swarthmore transit finder-chart tool](https://astro.swarthmore.edu/transits/finding_charts.cgi)

### Methods and advanced validation

- [Kreidberg 2015: `batman` transit model](https://arxiv.org/abs/1507.08285)
- [Kostov et al. 2019: DAVE vetting pipeline](https://arxiv.org/abs/1901.07459)
- [Giacalone et al. 2021: TRICERATOPS](https://arxiv.org/abs/2002.00691)
- [Gomez Barrientos et al. 2025: TRICERATOPS+ with ground-based light curves](https://arxiv.org/abs/2508.02782)
- [Morton et al. 2016: `vespa`](https://arxiv.org/abs/1605.02825) and [Morton & Johnson 2011: FPP framework](https://arxiv.org/abs/1101.5630)
- [Hadjigeorghiou et al. 2025: RAVEN methodology](https://arxiv.org/abs/2509.17645) and [Lafarga Magro et al. 2026: RAVEN TESS-SPOC FFI search](https://arxiv.org/abs/2603.22597) (tables: Zenodo 19661443) — cite only if TIC 390988385 membership is confirmed
- [Castro-González et al. 2024: the Neptunian desert/ridge/savanna](https://arxiv.org/abs/2409.10517)
- [TESS-cont contamination tool](https://github.com/castro-gzlz/TESS-cont)
- [ORBYTS student ephemeris-refinement study](https://arxiv.org/abs/2005.01684) and [TESS ephemeris-recovery study](https://arxiv.org/abs/1906.02197)
- Foreman-Mackey et al.: `emcee`; Speagle: `dynesty`; Foreman-Mackey et al.: `celerite2`; Hippke & Heller: TLS — cite the exact versions of whichever samplers/tools are actually adopted
- Young 1967 / Osborn et al. 2015: scintillation noise formulae
- SOAR TESS speckle survey (Ziegler et al. series) and ShARCS/SImMER documentation for the imaging provenance — identify the exact publication attached to the 2021-10-01 SOAR observation before citing

The final bibliography should cite the exact software versions, data releases, catalog queries, and scientific definitions actually used. This list is a starting set, not permission to cite a source without reading it.
