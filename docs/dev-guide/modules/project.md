# Projects Documentation

A New Project should be initiate for any kind of tools. This explains how does a New Project created along with its credentials and stored into database using `Project()` Module.

#####Technologies used by application:
Python version - 2.7

## Importing Modules

 Importing some of the external modules and sub-module in `rsquarelabs_core.config` and `rsquarelabs_core.engines.db_engine`
that are needed in this module.

```python
import os
from datetime import datetime
from time import time
from termcolor import cprint
from rsquarelabs_core.config import RSQ_PROJECTS_HOME, RSQ_DB_PATH
from rsquarelabs_core.engines.db_engine import DBEngine
```

Create the object for `DBEngine()`.
```python
db_object = DBEngine(RSQ_DB_PATH)
```

## How to Save a New Project

A class is created as `Project()` and arguments needed for this class are initiated in `__init__()` method. Here we sent
 arguments as dictionary.

```python
object = Project(project_title=project_data["title"], project_tags=project_data["tags"], project_user_email=project_data["user_email"],
                 project_short_note=project_data["short_note"], project_slug=project_data["slug"])
```

Project is saved into database using `object.save()`.

```python
def save(self):

    """
    This saves a new project details like log files into database using sqlite3 query. Some other methods are integrated
    into this method.
    :return: Returns the project identification number which is created.
    """
    db_object.do_insert("INSERT INTO projects (title, tags, user_email, slug, short_note, path, config, log, type, date, is_delete)\
            VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                                   (self.project_title,
                                    self.project_tags,
                                    self.project_user_email,
                                    self.project_slug,
                                    self.project_short_note,
                                    project_path,
                                    project_config,
                                    project_log,
                                    project_type,
                                    project_date,
                                    is_delete,))
            project_id = db_object.cur.lastrowid
            return project_id
```
#####Methods which are integrated into Save()

Creates the log file for the project.

```python
def create_log(self, project_log=None, project_data=None, project_id=None):
    """
    Creates the log file for the project.
    :param project_log: Path for the project log file.
    :param project_data: Data regarding on project for the log file.
    :param project_id: Project identification number for the log file.
    """
    fh_log = open(project_log, 'w', 0755)

```

Creates the configure file for the project.

```python
def create_config(self, project_config=None):
    """
    Creates the configure file for the project.
    :param project_config: Path for the project config file.
    """
    fh_config = open(project_config, 'w', 0755)
```

Generating the path for a project to be saved as project directory.

```python
def generate_path(self, slug=None):
    """
    This generates the path for a project to be saved as project directory.
    :param slug: Slug of a project for naming the project directory.
    :return: Returns project path that as generated.
    """

     if slug == None:
        slug = self.generate_slug()

    path = os.path.join(RSQ_PROJECTS_HOME, slug)

    os.mkdir(path, 0755)
    return path

```

Generating slug for a project.

```python
def generate_slug(self):
    """
    This generates slug for a project.
    :return:Returns slug.
    """

    slug = self.project_title.replace(" ","-").replace("_","-")\
            .replace("/","-").replace("\\","-").replace(".","-").replace(",","-").replace(";",'-')\
            .replace(":","-").replace("--","-")

    return slug
```

## How to create a New Run in a Project

Having runs to a project provides us good experience. The `object.save_run()` method saves the runs of a project into
database.
```python
object.save_run(run_name=run_name, version=version, run_data=run_data,
                              parent_run_id=parent_run_id, master_id=master_id,is_delete=is_delete,
                              protocol_class_name=protocol_class_name, project_id=project_id)
```
```python
def save_run(self, *args, **kwargs):
    """
    This saves run of a project into database.
    :param kwargs: credentials required for saving run into database.
    :return: Returns run's identification number of a project.
    """
    db_object.do_insert(" INSERT INTO runs (run_name, version, parent_run_id, master_id, run_data, is_delete, project_id, class_name)\
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)", (run_name, version, parent_run_id, master_id, run_data, is_delete, project_id, protocol_class_name, ))

    run_id = db_object.cur.lastrowid
    return run_id
```





