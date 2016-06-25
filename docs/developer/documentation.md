# Documentation




## Adding a Page

To add a new page for the documentation, please follow the steps: 

1. Create a page `some-page.md` in `docs/` ie., in the location `docs/some-page.md`. 

2. Add the page info to `mkdocs.yml` in the project root folder ie., `framework/mkdocs.yml`
```
- Dev Guides:
    - Environment Setup: developer/setting-environment.md
    - Some Page : developer/somepage.md

```
**Note:** Create the page in `developer` folder in `docs` if it belongs to developer docs, or in `user` for user docs

## Generating the Docs 

Use `mkdocs serve` , then you will be able to access the `http://localhost:8000` . 


## Deploying to gh-pages

Make the necessary changes and pull the latest code, then do the command 
```
mkdocs gh-deploy --clean
```
It will ask you for your github credentails, if you have permissions to push to the framework repository, it will generate a build and deploys it to the 
