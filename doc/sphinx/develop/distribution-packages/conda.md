# Distributing Conda Package via conda-forge

The Conda recipe hosted at <https://github.com/conda-forge/cantera-feedstock> is used
to build packages for both stable and pre-release versions of Cantera.

In the instructions below, `upstream` refers to the repository above, while `origin`
refers to your fork of this repository.

## Stable Cantera releases

- Fork the repository above and add it as a remote
- Add the remote repository used by the bot that updates dependencies. From the repo
  directory, run the commands:
  ```
  git remote add regro-cf-autotick-bot https://github.com/regro-cf-autotick-bot/cantera-feedstock
  git fetch --all
  ```
- Update `recipe/meta.yaml`:
  - Change the header line that looks like `{% set version = "3.1.0" %}`
  - In the `source` section, update the `url` and `sha256` for both the main `cantera`
    repo as well as the `cantera-example-data` repo.
    - The hash in the `url` is the hash of the commit corresponding to the tagged
      release.
    - The easiest way to get the `sha256` of the archive is to download it from that URL
      and run `sha256sum` on that file.
  - Update the value of `CT_GIT_COMMIT` with the same hash as used in the URL for the
    `cantera` repository.
  - Under `build`, reset the `number` to 0.
- Use `git cherry-pick` to apply changes from any PR branches opened by
  `regro-cf-autotick-bot`
- Do whatever else is needed to bring the conda-forge recipe up to date, for example
  running `conda-smithy rerender`.
- Create a pull request against the `main` branch of `conda-forge/cantera-feedstock`
- Once the PR is merged, you should be able to see the status of the final package
  builds in [Azure DevOps](https://dev.azure.com/conda-forge/feedstock-builds/_build?definitionId=11466&_a=summary)

## Cantera pre-releases

Whether your starting point should be the existing `main` or `dev` branch is a bit
tricky. In either case, there are likely to be conflicts to work out. Use caution to
make sure that important updates to the build aren't lost when merging branches.

### Option 1: Start from latest `main` branch

This process makes the most sense for the first `dev` build after a stable Cantera
release, or if the `main` branch has had significant updates since the last `dev` build.

- Check out the `dev` branch from `upstream`.
- Create a new working branch
- Run `git merge upstream/main` to bring in all changes from the `main` branch, likely
  preferring changes from `main` in case of merge conflicts.
- Follow the same steps as above for updating `recipe/meta.yaml`.
- Ensure that `recipe/conda_build_config.yaml` still contains the block:
  ```yaml
  channel_targets:
  - conda-forge cantera_dev
  ```
- Create the pull request against the `dev` branch of `conda-forge/cantera-feedstock`.

### Option 2: Start from the latest `dev` branch

This process is most useful when the `dev` branch is relatively up to date and there
aren't any relevant changes to the `main` branch.

- Check out the `dev` branch from `upstream`.
- Create a new working branch
- Follow the same steps as above for updating `recipe/meta.yaml`.
- Create the pull request against the `dev` branch of `conda-forge/cantera-feedstock`.
