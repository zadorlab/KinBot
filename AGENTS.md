# Working directives

- Keep KinBot general. Diagnose failures using native calculation evidence and
  attempt provenance, then implement a shared recovery procedure for the failure
  class. Check existing successful recovery strategies before proposing another
  one-off fix. Do not branch on example molecules, run paths, paper results, or
  production identifiers in installed code.
- Recovery procedures must check applicability, retain attempt history, avoid
  repeating exhausted settings, and use bounded retries. Keep accepted geometry,
  theory, electronic-state definitions, and scientific validation unchanged
  unless the user explicitly requests a refinement or model change.
- Use the owning drivers and existing live controls for targeted retries. Archive
  affected attempts and follow their dependency chains; preserve unrelated work.
  Distinguish implemented/tested recovery logic from a successful native retry.
- Keep paper comparisons and example-specific validation outside installed
  production code. Examples can provide regression evidence, not assumed answers.
- Before providing a controller-refresh handoff, test the complete orchestration
  as well as the changed helper: transient submission/staging races, preservation
  of live work, decision application, sampling adoption, and repeated invocation.
  Distinguish local simulated-site validation from a verified remote handoff.
- For a fixable run problem, carry diagnosis through implementation, validation,
  publication, and concrete pull/update/refresh commands. Diagnosis alone is not
  a completed task. If evidence or access prevents a fix, state that specific
  blocker rather than presenting an unverified change as a solution.
