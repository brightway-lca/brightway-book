# Brightway Documentation (JupyterBooks)

[![Join the chat at https://gitter.im/brightway-lca/documentation](https://badges.gitter.im/brightway-lca/documentation.svg)](https://gitter.im/brightway-lca/documentation?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/brightway-lca/brightway-documentation-jupyter-books/discussions)

## ⁉️ FAQ

_"Where can I find the developer documentation?"_ \
The API documentation is at [`brightway.readthedocs.org`](https://brightway.readthedocs.org/).

_"How can I contribute to the documentation?"_ \
If you are familiar with Jupyter Book, you can add your changes to a fork of this repo and open a pull request. \
If you would rather get help from a developer, please [start a new discussion here](https://github.com/brightway-lca/brightway-documentation-readthedocs/discussions). Your suggestions/additions will be added for you.

_"The live code execution in the Jupyter Book is broken!"_ \
C'est La Vie

_"Where can I get help beyond the documentation?"_ \
Use the Gitter channel linked above.

## ☑️ Instructions

1. clone the repository
2. follow the [`jupyter-book` setup instructions](https://jupyterbook.org/en/stable/start/overview.html)
3. follow the [`ghp-import` setup instructions](https://jupyterbook.org/en/stable/start/publish.html)
4. edit content or add new files
5. [build the book](https://jupyterbook.org/en/stable/start/build.html) by running, from within the repository folder

```
jupyter-book build ./
```

6. preview the book by opening

```
_build/html/index.html
```

7. [publish the book](https://jupyterbook.org/en/stable/start/publish.html) by running, from within the repository folder

```
ghp-import -n -p -f _build/html
```

8. check if the [`gh-pages` branch](https://github.com/brightway-lca/brightway-documentation-jupyter-book/tree/gh-pages) of the repository has been updated
9. check if the [`documentation.brightway.dev`](https://documentation.brightway.dev/) website has been updated

## 📚 References

Compare the `jupyter-book`:

1. [documentation](https://jupyterbook.org/en/stable/intro.html)
2. [feature requests queue](https://executablebooks.org/en/latest/feature-vote.html)
3. [discussions on GitHub](https://github.com/orgs/executablebooks/discussions)