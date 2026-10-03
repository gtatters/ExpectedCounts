# Expected Counts in Contingency Tables

An interactive Shiny app for BIOL 3P96 (Biostatistics) at Brock University.

## What this app does

When we analyse a contingency table (for example, with a chi-square test), we
compare the counts we observed with the counts we would *expect* if the two
traits had nothing to do with each other. A common first guess is that the
expected count in each box should simply be the total divided by the number of
boxes (N / 4 for a 2 × 2 table). This app shows why that guess is wrong, and
where the correct expected counts (Row total × Column total / N) come from.

The app uses a single hypothetical example throughout: frogs sampled from
forest and meadow ponds, each scored as infected or not infected with a skin
fungus. All data are simulated for teaching purposes.

Three tabs build the idea step by step:

1. **Where expected counts come from** — A frequency tree splits the frogs by
   habitat and then by infection status, showing that each expected count is
   N × P(row) × P(column). This is the Multiplication Rule for independent
   events. A mosaic plot shows the same expected counts as areas. An optional
   overlay shows the N / 4 guess, which only matches when both traits are split
   exactly 50:50.
2. **One sample** — You choose the true process in nature (no link between
   habitat and infection, or a real link). The app then samples frogs at random
   and shows the observed table (one dot per frog), the expected table built
   from R × C / N, and the N / 4 guess side by side. The chi-square statistic is
   calculated both ways for comparison.
3. **Many samples** — The study is repeated hundreds of times. A histogram
   shows that observed counts average out to R × C / N, not N / 4, when there is
   no link. A bar chart shows how often a chi-square test claims a link using
   each set of expected counts.

## Key ideas illustrated

- Expected count = N × P(row) × P(column) = Row total × Column total / N
- Expected counts keep the same row and column totals as the observed data;
  they only redistribute the same individuals among the boxes
- "No link" looks like a flat dividing line across columns in a mosaic plot
- Using N / 4 answers a different question ("are all four combinations equally
  common?"), so it raises false alarms whenever the groups differ in size

## How to use

1. Start on Tab 1 and move the habitat and infection sliders. Tick the
   **equal split** box to compare the correct expected counts with N / 4.
2. On Tab 2, choose **No link** or **Linked**, set the sample size and
   percentages, and press **Resample** to catch a new set of frogs.
3. On Tab 3, use the same settings and press **Resample** to repeat the whole
   study many times. Compare how often each method claims a link.

## Learning goals

- Explain why expected counts are not simply N divided by the number of boxes
- Connect the R × C / N formula to the Multiplication Rule for independent events
- Recognise what independence looks like in a table and in a mosaic plot
- Understand why using the wrong expected counts leads to false conclusions
  in a chi-square test

## Notes

- Chi-square statistics in the app are calculated without a continuity
  correction, so values may differ slightly from R's `chisq.test()`, which
  applies Yates' correction to 2 × 2 tables by default.
- Companion app: **Probability: Concepts and Rules** (see the Multiplication
  Rule tab).

## Course context

Developed for BIOL 3P96 — Biostatistics, Brock University.
Built with R and Shiny (base R graphics only).
