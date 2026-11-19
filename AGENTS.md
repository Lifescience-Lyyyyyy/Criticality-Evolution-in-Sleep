# Project instructions

- Use GitHub Actions MATLAB R2024b for the pinned Sooter estimator. Inspect existing workflows before editing them. Run one stage at a time and require a passing scientific gate before starting dependent experiments.
- Run local development and tests through WSL. Use Colab for data preparation and Drive for large inputs, intermediate arrays, and durable archives. Keep credentials, raw recordings, and large binaries out of Git and Actions artifacts.
- Validate external inputs and prerequisites once at the entry point. Rely on established invariants downstream, surface unexpected errors with diagnostics, and add retries or recovery for known failure modes. Enforce scientific QC and protect data integrity.
- Preserve upstream algorithms and published analysis choices. Record code revisions, source file checksums, experiment parameters, results, and documented deviations. Use the animal or recording as the appropriate experimental unit.
