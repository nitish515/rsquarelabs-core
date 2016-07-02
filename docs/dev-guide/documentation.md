# Documentation Rules

This section explains how to document the modules inside this tool.


## Adding a Page

To add a new page for the documentation, please follow the steps: 

1. Create a page `some-page.md` in `docs/` ie., in the location `docs/some-page.md`. 
Add in an appropriate section ie., `user-guide` or `dev-guide`.

2. Add the page info to `mkdocs.yml` in the project root folder ie., `core-client/mkdocs.yml`
```
- Dev Guides:
    - Environment Setup: developer/setting-environment.md
    - Some Page : dev-guide/somepage.md

```
**Note:** Create the page in `dev-guide` folder in `docs` if it belongs to developer docs, or in `user-guide` for user docs

## Generating the Docs 

Use `mkdocs serve` , then you will be able to access the `http://localhost:8000` . 

You can use `mkdocs build` to generate a build in the folder `site` .


## Deploying to gh-pages

Make the necessary changes and pull the latest code, then do the command 
```
mkdocs gh-deploy --clean
```
This will create a build in the folder `site` and push to the `gh-pages` branch.

It will ask you for your github credentails, if you have permissions to 
push to the core-client repository, it will generate a build and deploys 
it to the gh-pages branch
