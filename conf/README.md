## Configuration Changes: `base_._original.config` vs. `base.config`

This section details the differences in resource allocation between the original configuration (`base_.original.config`) and the modified version (`base.config`). The primary changes involve reductions in CPU allocations for specific process labels.

**Summary of Changes:**

| Process Label         | Setting | Original Value (`base.copy.config`) | Modified Value (`base.config`)   | Change    |
| :-------------------- | :------ | :------------------------------------ | :------------------------------- | :-------- |
| `process_low`         | `cpus`  | `{ check_max( 8 * task.attempt, 'cpus' ) }` | `{ check_max( 4 * task.attempt, 'cpus' ) }` | Decreased |
| `process_medium`      | `cpus`  | `{ check_max( 8 * task.attempt, 'cpus' ) }` | `{ check_max( 4 * task.attempt, 'cpus' ) }` | Decreased |
| `process_high_disk`   | `cpus`  | `{ check_max( 16 * task.attempt, 'cpus' ) }` | `{ check_max( 4 * task.attempt, 'cpus' ) }` | Decreased |
| `filter_reads`        | `cpus`  | `8`                                   | `4`                              | Decreased |
| `mapBothReads`        | `cpus`  | `{ check_max( 16 * task.attempt, 'cpus' ) }` | `{ check_max( 4 * task.attempt, 'cpus' ) }` | Decreased |

**Key Observations:**

* CPU resources have been reduced for processes labeled `process_low`, `process_medium`, `process_high_disk`, `filter_reads`, and `mapBothReads` in the `base.config` file.
* Other settings, including default resources (memory, time), error handling strategies, and configurations for other labels (`process_high`, `process_long`, `process_high_memory`, `extract_reads`, `fastqc`), remain unchanged between the two files.

This adjustment likely aims to optimize resource usage based on observed requirements or cluster availability.