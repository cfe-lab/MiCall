## BaseSpace

This section explains how to publish and test MiCall as an **Illumina
BaseSpace Native App**.

**BaseSpace** is Illumina's platform for running analysis applications
against sequencing data. A **Native App** in BaseSpace is a
containerized application that BaseSpace can launch for a user. For
MiCall, the important BaseSpace-side pieces are:

* **App**: the BaseSpace app entry that users see and launch.
* **Form**: the configuration page shown to the user before
  launch. This is where inputs, options, and launch behavior are
  defined.
* **Form Builder**: the BaseSpace web UI where developers edit the
  form and its launch configuration. Some Illumina documentation calls
  this **Forms Builder**.
* **Report**: the BaseSpace-side presentation layer for results, if
  one is configured.
* ``: the launch specification used by BaseSpace to decide which
  container image to run, what command to run inside it, and which
  runtime options to enable.
* **SpaceDock**: the BaseSpace agent that receives jobs, stages data,
  launches the app container, and uploads results back to BaseSpace.
* **Local agent**: a developer-run SpaceDock instance used for testing
  with **Send to Local Agent**.
* **Continue**: the normal cloud execution path inside Illumina's
  infrastructure.

MiCall's BaseSpace entrypoint is the `basespace` subcommand. The
launch spec in `callbacks.js` tells BaseSpace to run that subcommand
inside the MiCall container.

This section assumes you are **publishing a new MiCall release**, not
setting up the whole MiCall development environment from scratch.

### Accounts and access you need

You will need:

* the shared Illumina account used for MiCall BaseSpace work,
* permission under that account to edit the MiCall app in BaseSpace,
* permission to push to `docker.illumina.com`,
* a built MiCall Docker image for the new release.

### Build and push the MiCall image

Build the MiCall Docker image on your normal host machine or in
CI. Push it to `docker.illumina.com` and record the resulting
digest. Pin the BaseSpace Native App to a specific Docker digest
rather than using a mutable tag such as `:latest`.

Log in to the private registry on the machine that will build and push
the image:

```shell
docker login docker.illumina.com
```

After pushing, record the digest:

```shell
docker images --digests=true | grep micall
```

You will use that digest in the BaseSpace launch configuration.

### How to get to the MiCall app in BaseSpace

Open the BaseSpace developer website in your browser and sign in with
the shared MiCall account.

Developer website:

```text
https://developer.basespace.illumina.com
```

Then navigate as follows:

1. Click **My Apps**.
2. Find and open the MiCall application.
3. Use the tabs under the app name to move between **Details**, **Form
   Builder**, and **Reports Builder**.

Most of the BaseSpace release work for MiCall happens inside that one
app page.

### Update the launch configuration in BaseSpace

From the MiCall app page in the Developer Portal:

1. Open the **Form Builder** tab.
2. In the revisions list, select the currently active editable form
   revision, or duplicate the current revision if you want a fresh
   revision for the release.
3. Open that revision.
4. In the editor, use the **Current Template** dropdown and select
   `callbacks.js`.

The `callbacks.js` template tells BaseSpace how to launch MiCall.

For MiCall, the important launch behavior is:

* command line: `['basespace']`
* BSFS enabled: `Options: ['bsfs.enabled=true']`
* container image: the new MiCall image, pinned by digest

Keep BSFS enabled unless there is a specific reason to test the
deprecated pre-download path.

Example:

```javascript
function launchSpec(dataProvider) {
    return {
        commandLine: ["basespace"],
        containerImageId: "docker.illumina.com/cfe_lab/micall@sha256:REPLACE_ME",
        Options: ["bsfs.enabled=true"]
    };
}
```

After updating the form revision:

1. Save the revision.
2. Return to the MiCall app page.
3. Create a new BaseSpace app version for the new MiCall release
   number.

If MiCall also has a report revision that should change for the
release:

1. Open the **Reports Builder** tab for the same app.
2. Open the relevant report revision there.
3. Update it.
4. Activate it when ready.

### Test MiCall through BaseSpace

The primary release-validation path should be **cloud execution**, not
a local agent. For the normal release flow, use **Continue** in
BaseSpace so that MiCall runs in Illumina's own execution
environment. This avoids installing SpaceDock or BSFS on your
development machine and tests the path that real users actually
receive.

To launch a cloud test run:

1. Open the MiCall app page in the Developer Portal.
2. Open the **Form Builder** tab.
3. Open the new form revision.
4. Fill in the form inputs for a test run and choose the output
   project.
5. Click **Continue**.
6. Wait for the BaseSpace cloud run to complete.
7. Review logs, outputs, and result rendering in BaseSpace.
8. Process a small set of MiCall microtests this way.

This is the main end-to-end test for the BaseSpace launch path. It
validates that:

* the BaseSpace form is wired correctly,
* `callbacks.js` points at the correct MiCall image,
* the `basespace` entrypoint works,
* BaseSpace can stage inputs and collect outputs,
* the results are accepted by BaseSpace,
* the cloud execution path behaves correctly for real users.

### Validate the cloud path carefully

Because **Continue** is the primary release-validation path here, take
the cloud run seriously. Check not only whether the run succeeds, but
also whether the outputs, logs, result layout, and user-visible
behavior match expectations.

This is especially important for catching differences in staged data,
environment behavior, permissions, available reference data, and
report rendering.

### Optional: isolated local-agent debugging

A local SpaceDock agent is optional. Use it only when you specifically
need to debug staging, local execution behavior, or
BaseSpace-to-container handoff in more detail than the cloud run
exposes.

Do **not** treat this as something that should live on your normal
development machine. If local-agent testing is needed, use a dedicated
Linux VM or separate Linux machine.

On that isolated machine:

* install Docker,
* install Mono,
* install the `spacedock` package,
* install the `bsfs` package,
* create `/data` and `/genomes`.

To start a local-agent session from BaseSpace:

1. Open the MiCall app page in the Developer Portal.
2. Open the **Form Builder** tab.
3. Open the form revision you want to test.
4. On the right side of the editor, find the sample command shown for
   starting a local agent.
5. Copy that command.
6. Run it on the isolated Linux machine.
7. Return to the Form Builder and use **Send to Local Agent** if you
   want to run the job through that machine.

The command will look roughly like this:

```shell
sudo spacedock -a <agent_id> -m <mission_control_uri>
```

If startup is successful, SpaceDock should report that it has
connected to the Docker service and has begun polling BaseSpace for
jobs for that application.

#### Why this is optional and isolated

* The Dockerized `basespace/spacedock` image is not a reliable default
  on modern Docker because it is very old and may fail to pull due to
  legacy image-manifest format issues.
* Installing SpaceDock and BSFS directly on your day-to-day
  development machine is often undesirable.
* **Continue** already tests the real cloud execution path, which is
  the most important validation path for release readiness.

If you do use local-agent testing, treat it as a debugging supplement
rather than the main release gate.

### Activate the release in BaseSpace

Once the new release has passed cloud-path testing, and optional
local-agent debugging if you chose to do it:

1. Open the MiCall app page in the Developer Portal.
2. Open the **Form Builder** tab and activate the new form revision.
3. If a report changed, open the **Reports Builder** tab and activate
   the new report revision.
4. Confirm that the new app version points at the intended revision
   set.
5. Verify that the user-visible version number matches the MiCall
   release number.

After that, follow the appropriate sharing or publication path for the
app's ownership and distribution model.

### Clean up after testing

If you used a separate local-agent machine, local runs will leave data
on that host. Clean old data from `/data` and remove unnecessary
Docker images and containers after testing. Keep `/data` around only
when you want to reuse staged inputs for debugging.

### Practical notes and gotchas

* Keep the currently active BaseSpace revision untouched until the new
  one has passed testing.
* Treat `callbacks.js` as release-critical: a wrong digest or wrong
  command line can make an otherwise good MiCall image fail in
  BaseSpace.
* If image pulls fail in BaseSpace-related testing, log out and log
  back in to `docker.illumina.com` on the machine that is pushing the
  image, then verify that the pushed digest is the one referenced in
  `callbacks.js`.
* If cloud execution behaves differently than expected, inspect
  whether the issue is in the MiCall image itself, the BaseSpace form
  wiring, or the result or report configuration.
* If you use local-agent debugging, remember that it is supplementary:
  differences between local staging and cloud staging are possible,
  especially around `/genomes`, permissions, and environment details.
* Do not depend on the `basespace/spacedock` Docker Hub image as the
  primary path on a modern system.
