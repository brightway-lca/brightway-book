# Brightway Interactive Learning (based on [Jupyter Book](https://jupyterbook.org/en/stable/intro.html))

[![Brightway](https://img.shields.io/static/v1?label=Brightway&message=ecosystem&color=45bfb0&logo=data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMTA2NSIgaGVpZ2h0PSI2OTAiIHZlcnNpb249IjEuMSIgdmlld0JveD0iMCAwIDEwNjUgNjkwIiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciPgogPGRlZnM+CiAgPGNsaXBQYXRoIGlkPSJjbGlwUGF0aDIxNzMiPgogICA8cGF0aCBkPSJtLTU5NSA0NDBoMWUzdi0xZTNoLTFlM3oiLz4KICA8L2NsaXBQYXRoPgogPC9kZWZzPgogPGcgdHJhbnNmb3JtPSJtYXRyaXgoMS4zIDAgMCAtMS4zIDY1MyA0MDMpIiBjbGlwLXBhdGg9InVybCgjY2xpcFBhdGgyMTczKSIgZmlsbD0iI2ZmZiI+CiAgPHBhdGggZD0ibTAgMGMwLTEuMi0xLjUtMy40LTAuNDctNS4yIDEuMy0yLjQgNS40LTEuNSA1LjgtM3YtMC4wMThsLTQuNy0wLjA1NC0yLTQuMWMtMS45IDAuNzEtMi40IDMuMi0zLjIgNi42LTAuNDQgMi0wLjcxIDMuNCAwLjA4OSA0LjctMTQtMy4zLTMwLTUuNC00NS01LjQtNDkgMC4wMzYtNzEgMTctMTA3IDI0IDEuNS0xLjEgMi43LTMuMyAyLjEtNC45bC02LjEgMi4yLTQtMy4xYy0xLjQgMS40LTAuNTggMi44LTIuMiA0LjgtMi4xIDIuNS01LjMgMi42LTUuMiAzIDAuODctMC4wNzIgMS43LTAuMTQgMi42LTAuMjItMS42IDAuMjQtMi41IDAuNC0yLjYgMC4yMi0zLjkgMC4zMS04IDAuNDktMTIgMC41MS02NiAwLjIyLTEwMi00Mi0xMjItNzkgMy42IDIuOSA4LjcgMS41IDEwIDEuMSA0LjItMS4xIDguOC00LjEgOC40LTYuMWwtNC41LTEuNSAwLjUyLTAuNTRjLTAuMjItMC4wMTktMC40My0wLjAxOS0wLjY1LTAuMDE5LTQuMiAwLjA1NS01LjEgMy45LTkuMiA0LjItMy44IDAuMjktNi40LTIuNy03LjItMS45LTEzLTI3LTE4LTUwLTE4LTUwaC0zNHM2LjcgMzQgMzAgNjhjLTEuNiAxIDEgNS0xLjEgOC4zLTIuNSAzLjgtOC40IDIuNC04LjggNC41LTAuMDE4IDAuMTMtMC4wMzYgMC4yNSAwIDAuMzhsNy4yLTAuNTIgMi41IDUuNWMxLjEtMC40IDItMS43IDMuOS00LjIgMS40LTEuOCAyLjQtMy4yIDMuMS00LjMgMjYgMzQgNjkgNjcgMTM5IDY3IDIuNSAwIDUtMC4wNzMgNy40LTAuMi0wLjU4IDIuMSA1IDMuNSA1LjcgOC4xIDAuNzQgNC44LTQuNyA3LjgtMy40IDEwIDAuMDE4IDAuMDE4IDAuMDM2IDAuMDU0IDAuMDU0IDAuMDcybDUuOC00LjQgNC43IDMuMWMwLjU2LTAuODcgMC42My0xLjctMC40NS04LTEuNC04LjItMS45LTktMi43LTkuNyA0MS01LjIgNjgtMjggMTE4LTI5IDE1LTAuNTQgMzAgMC43OCA0NCAzLjEtMC4yNyAwLjE2LTAuNTIgMC4zNi0wLjc2IDAuNjItMSAxLjEtMS4xIDMuMS0xLjQgNy4xLTAuMjUgNC4zLTAuNCA2LjYgMC40MiA3LjdsMi44LTMuMiAzLjIgMS4yYzAuMDM2LTAuMDcyIDAuMDU0LTAuMTYgMC4wNzItMC4yNCAwLjQ0LTEuOC0xLjEtMi43LTEuMS01LjYgMC0zLjUgMi4zLTQuOSAxLjgtNi41LTAuMDM2LTAuMDktMC4wNzItMC4xOC0wLjExLTAuMjUgNDUgOC43IDc4IDI4IDc5IDI4IDAtMC42My0zNS0yMi03OS0zM20tMzMyLTMwLTU1IDQuMS0zNSA0MC0xNyAyOC0xNSAyNC05LjcgMjctMjYgMTYgMzgtOS40IDM5LTI0IDMxLTI4IDI3LTM5IDE5LTMzem00MTEgNjUgMTEgNzkgNTcgNDYgNjggOS4zIDQ1LTAuMzkgNDcgMTYtMTYtMzYtMTMtMzQtMjQtNTgtMzItNDktMzktMjAtNDMgMy43em04OC0yNDktMTggMzItMjEgMjYtMzUgMzQtMjEgMjYtMjAgMjItMzcgNDgtMTcgMjAgNzAgMC44NyA1Ni00OSAyMi00MCAxNS0zNyAyLjgtMzJ6bTAgMC0xNyAxNy0zOSA2LjItNDQgMTMtNTAgMzMtMzAgNTUtMy4zIDUzIDEzIDI3IDEuOSAzLjcgMTctMjAgMzctNDggMjAtMjIgMjEtMjYgMzUtMzQgMjEtMjZ6bS05MSAzNDgtMTYtNTctMzEtNDYtMzMtMTIgMC4xNSAzLjkgMS42IDQ0IDQuNCAzNiA1LjEgMzQgNC4zIDQwIDIgMzQgMi4zIDM2IDEzIDQyIDQuNS0zNyAxNi0yMiAyMC0zN3ptLTQ3IDE1NC0xMy00Mi0yLjMtMzYtMi0zNC00LjMtNDAtNS4xLTM0LTQuNC0zNi0xLjYtNDQtMC4xNS0zLjktMi44IDMuMi00OCA1NS03LjIgMzcgMC43OCA2NCAyNiA0OCAzMyAyOXptLTE0NS0zOTktMjAgNDMtMTIgNDEtNi4xIDI0LTYuNyAyMCA1OC0yMSAxNS0yNiAxNi00NC00LjItMzctMTItNjEtNS4yIDI0em0yOC02MS0yNyAxOS0yMCAxMS0yOCAyMC0zMSAzMC02LjggNDggMTcgNDIgMjQgMTggMS41LTQuNiA1LjItMTYgNi4xLTI0IDEyLTQxIDIwLTQzIDIyLTM3em0tMTEgMzA1LTE4LTUxLTUyLTM0LTAuMjUgNS44LTAuODMgMTkgMS4yIDM5LTMuMSA0My02LjcgMzgtMTEgNDYtMTUgNjMgMjAtMjUgNDMtNDEgMzAtNDl6bS03MC03OSAwLjI1LTUuOC01NiA0Mi0xNyA2NSAyLjcgNTAgOS40IDQwIDE2IDE3IDguNCA0MCAxNS02MyAxMS00NiA2LjctMzggMy4xLTQzLTEuMi0zOXptLTU0LTE5NC0xMiAyMC0yMiAyNi0xOCAxNS0xNyAxNy0wLjEzLTAuNTYtOS43LTQ1IDIxLTM0IDIzLTEzIDIxLTguMiAzNy0xOS0xNSAxM3ptMjQtNDEtMTUgMTMtOC4zIDI4LTEyIDIwLTIyIDI2LTE4IDE1LTE3IDE3IDQyIDE0IDI0LTQuOCAxNS0xNiA4LjYtMzMgMS41LTI4LTMuNC0yMHptLTExMCAyMDItMjMtNTEtMi44IDQuOS0xOSAzMy0yNyAzOS0zMSAyOC0zOSAyNC0zOCA5LjQgMzIgMS4zIDMxIDEuNiAzMC0wLjI5IDQzLTE2IDMwLTMweiIvPgogPC9nPgo8L3N2Zz4K)](https://github.com/brightway-lca)
![License](https://img.shields.io/github/license/brightway-lca/brightway-learn?color=green&logo=Open%20Source%20Initiative&logoColor=white) [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat&logo=GitHub)](https://github.com/brightway-lca/brightway-documentation-readthedocs/discussions) \
Maintainance Team: [@michaelweinold](https://github.com/michaelweinold)

The interactive learning and teaching documentation for the Brightway life-cycle assessment software package. [Readthedocs.org](https://about.readthedocs.com/?ref=readthedocs.com) with the [Jupyter Book Theme](https://jupyterbook.org/en/stable/intro.html) and the [Thebe tool](https://thebe.readthedocs.io/en/stable/) are used to build/host the interactive documentation.

| static documentation | interactive documentation | development playground |
| ---- | ------------- | ------------- |
| [docs.brightway.dev](https://github.com/brightway-lca/brightway-documentation) | [learn.brightway.dev](https://github.com/brightway-lca/brightway-book) | [live.brightway.dev](https://github.com/brightway-lca/brightway-live) | 

## Quickstart

### Setup Repository

1. [Clone this repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository):

```bash
git clone https://github.com/brightway-lca/brightway-book.git
```

### Setup Python Environment

Set up a Python virtual environment that includes all packages required to build the documentation. A [Conda environment file](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) is provided [for convenient setup](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file). The file is located at [``./environment.yml``](environment.yml). Install the environment `jbook` by running from the repository root directory:

```bash
conda env create -f environment.yml
```

and activate the environment:

```bash
conda activate jbook
```

### Building the Documentation

> [!TIP]
> You can check [the Jupyter Book documentation](https://jupyterbook.org/en/stable/intro.html) for further information on how to structure your content.

1. You can build the documentation [by running the `build` command](https://jupyterbook.org/en/stable/start/build.html):

```bash
jupyter-book build ./
```

2. Preview the built documentation by opening the generated `html` file:

```
_build/html/index.html
```

### Publishing the Documentation

1. [Create a pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) to have a maintainer merge your changes into the `main` branch. As soon aschanges have been merged, Readthedocs.org will build the documentation and deploy it.

## 📚 References

Compare the `jupyter-book`:

1. [documentation](https://jupyterbook.org/en/stable/intro.html)
2. [feature requests queue](https://executablebooks.org/en/latest/feature-vote.html)
3. [discussions on GitHub](https://github.com/orgs/executablebooks/discussions)
