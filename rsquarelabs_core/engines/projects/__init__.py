from rsquarelabs_core.config import RSQ_PROJECTS_HOME, RSQ_DB_PATH
import os
from rsquarelabs_core.engines.db_engine import DBEngine
from datetime import datetime


db_object = DBEngine(RSQ_DB_PATH)

class Project(object):

    def __init__(self):
        pass

    def create(self, *args, **kwargs):
        # check if the project-key exist in db,
        # check if the project-key exist in rsquareProjects

        self.project_title = kwargs.get('project_title', None)
        self.project_tags = kwargs.get('project_tags', None)
        self.project_user_email = kwargs.get('project_user_email', None)
        self.project_slug = kwargs.get('project_slug', None)
        self.project_short_note = kwargs.get('project_short_note', None)

        PROJECT_PATH = os.path.join(RSQ_PROJECTS_HOME, self.project_slug)

        if os.path.exists(PROJECT_PATH):
            print "ERROR: Project name already exists"
            os.rmdir(PROJECT_PATH)
            os.mkdir(PROJECT_PATH, 0755)
        else:
            os.mkdir(PROJECT_PATH, 0755)

        project_type ="r2_gromacs"
        project_date = datetime.now().strftime("%Y-%m-%d %H:%M")
        is_delete = 0
        project_path = PROJECT_PATH
        project_log = os.path.join(PROJECT_PATH, 'r2_gromacs.log')
        project_config = os.path.join(PROJECT_PATH, 'r2_gromacs.config')

        cur = db_object.do_insert("INSERT INTO projects (title, tags, user_email, slug, short_note, path, config, log, type, date, is_delete)\
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
        return cur