# Project Overview

This repository contains all materials related to two radar signal processing quizzes completed as part of the course "Radar Systems". It includes code, reports, processed results, and inferences for both Quiz 3 and Quiz 4.

---

## Repository Structure

```
C:\Users\aarya\Github_Projects\Radar_Systems_Matlab_Quizzes
│   README.md
│
├───2022006_RS-Quiz-3
│       2022006_RS_Quiz-3.pdf
│       main.m
│       matched_filter.png
│       received_rx.png
│       transmitted_tx.png
└───2022006_RS-Quiz-4
        Aarya_Gupta_2022006.m
        Aarya_Gupta_2022006_HTML_Published.pdf
        Aarya_Gupta_2022006_method2.m
        Aarya_Gupta_2022006_PDF_Published.pdf
        Aarya_Gupta_2022006_Theory.pdf
        method2_report.pdf
        rs quiz 4.pdf
```

---

## Quizzes Overview

### Quiz 3

* **Objective:** Simulate and process automotive radar signal using a polyphase code to generate range-Doppler maps.
* **Contents:**

  * MATLAB code defining radar parameters, target properties, matched filtering, and Doppler processing.
  * Plots illustrating transmit pulse, matched filter response, and range-Doppler visualization.
  * Written answers to theory questions addressing resolution, ambiguity, and target detectability.

### Quiz 4

* **Objective:** Perform a similar simulation for a coded automotive radar (polyphase sequence) and analyze target resolution in range and Doppler domains.
* **Contents:**

  * MATLAB script for generating baseband signal, simulating received echoes for two targets, applying pulse compression, and computing the range-Doppler map.
  * Code comments and written calculations for parameters such as range resolution, Doppler resolution, and unambiguous limits.
  * Comparative discussion on whether the targets are resolved in range, Doppler, and both, as well as their placement within unambiguous radar limits.

---

## How to Use

1. **MATLAB Environment:** Ensure you have MATLAB R2024b (or later) installed.
2. **Run Simulation:** Navigate to the `src/` directory and open `quiz3_simulation.m` or `quiz4_simulation.m`. Execute the script.
3. **View Results:** Generated figures will be saved in the `results/` folder. Open the PDF reports for detailed discussions and inference.
4. **Modify Parameters:** To test different KEY values or target scenarios, edit the `KEY`, `R1`, `v1`, `R2`, and `v2` variables in the scripts.

---

## Inferences & Findings

* **Range Resolution:** Verified theoretical and simulated resolutions match (≈ 5.625 m for 16-bit polyphase code).
* **Doppler Resolution:** Calculated Doppler resolution (\~2.54 m/s) and confirmed via FFT processing.
* **Target Resolution:** Both quizzes demonstrate that the two targets at 100 m and 50 m, with velocities of ±5 m/s and ±2.5 m/s, are resolvable in range and Doppler given the selected parameters.
* **Unambiguous Limits:** Both targets fall well within maximum unambiguous range (900 m) and velocity (±162.34 m/s).

---

## Contact

For questions or feedback, please reach out to Aarya Gupta at \[[your.email@example.com](mailto:aarya22006@iiitd.ac.in)].
