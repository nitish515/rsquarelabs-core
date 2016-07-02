import os
from datetime import datetime
from time import time
from termcolor import cprint
from rsquarelabs_core.config import RSQ_PROJECTS_HOME, RSQ_DB_PATH
from rsquarelabs_core.engines.db_engine import DBEngine


__VERSION__ = "0.1dev"

db_object = DBEngine(RSQ_DB_PATH)

class Project(object):

    def __init__(self, *args, **kwargs):
        self.project_title = kwargs.get('project_title', None)
        self.project_tags = kwargs.get('project_tags', None)
        self.project_user_email = kwargs.get('project_user_email', None)
        self.project_short_note = kwargs.get('project_short_note', None)

        self.project_slug = kwargs.get('project_slug', None)
        if self.project_slug is None and self.project_title is not None:
            self.project_slug = self.generate_slug()





    def save(self):
        """
        This saves a new project details like log files into database using sqlite3 query. Some other methods are integrated
        into this method.
        :return: Returns the project identification number which is created.
        """

        project_path = self.generate_path(slug=self.project_slug)

        project_type ="r2_gromacs"
        project_date = datetime.now().strftime("%Y-%m-%d %H:%M")
        project_log = os.path.join(project_path, 'r2_gromacs.log')
        project_config = os.path.join(project_path, 'r2_gromacs.config')

        project_data = {}
        project_data["title"] = self.project_title
        project_data["tags"] = self.project_tags
        project_data["user_email"] = self.project_user_email
        project_data["short_note"] = self.project_short_note
        project_data["slug"] = self.project_slug
        project_data["path"] = project_path
        project_data["type"] = project_type
        project_data["date"] = project_date
        project_data["log"] = project_log
        project_data["config"] = project_config
        is_delete = 0


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

        # return cur
        project_id = db_object.cur.lastrowid

        # print project_id
        self.create_log(project_log=project_log, project_data=project_data, project_id=project_id)
        self.create_config(project_config=project_config)

        return project_id

    def create_log(self, project_log=None, project_data=None, project_id=None):
        """
        Creates the log file for the project.
        :param project_log: Path for the project log file.
        :param project_data: Data regarding on project for the log file.
        :param project_id: Project identification number for the log file.
        :return:
        """
        # print "----->4"
        fh_log = open(project_log, 'w', 0755)
        from random import randint
        project_create_details = project_data  # json.loads(project_data)
        project_create_details['project_id'] = randint(1, 1000)
        fh_log.write(
            "# RSQUARELABS-CORE v%s \n# Written by Ravi RT Merugu \n# https://github.com/rsquarelabs/rsquarelabs-core\n\n\n" % __VERSION__)

        mesg = """                  ============================================
                        Project created with id '%s',
                  ============================================""" % project_id
        # fh_config.write(cur.lastrowid)
        cprint(mesg, "green")

    def create_config(self, project_config=None):
        """
        Creates the configure file for the project.
        :param project_config: Path for the project config file.
        """
        # print "------>5"
        fh_config = open(project_config, 'w', 0755)



    def save_run(self, *args, **kwargs):
        """
        This saves run of a project into database.
        :param kwargs: credentials required for saving run into database.
        :return: Returns run's identification number of a project.
        """
        run_name = kwargs.get('run_name', "Default")
        version = kwargs.get('version', '1')
        parent_run_id = kwargs.get('parent_run_id', '0')
        master_id = kwargs.get('master_id', '0')
        run_data = kwargs.get('run_data', "Default")
        is_delete = kwargs.get('is_delete', '0')
        project_id = kwargs.get('project_id', None)
        protocol_class_name = kwargs.get('protocol_class_name', "Default")

        db_object.do_insert(" INSERT INTO runs (run_name, version, parent_run_id, master_id, run_data, is_delete, project_id, class_name)\
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?)", (run_name, version, parent_run_id, master_id, run_data, is_delete, project_id, protocol_class_name, ))

        run_id = db_object.cur.lastrowid

        return run_id


    def generate_path(self, slug=None):
        """
        This generates the path for a project to be saved as project directory.
        :param slug: Slug of a project for naming the project directory.
        :return: Returns project path that as generated.
        """
        # print "------>2"
        if slug == None:
            slug = self.generate_slug()

        path = os.path.join(RSQ_PROJECTS_HOME, slug)

        if os.path.exists(path):
            backup_path = path + "_%s" % time()
            # backup_slug = slug + "_%s" % time()
            os.rename(path, backup_path)


        os.mkdir(path, 0755)
        return path

    def generate_slug(self):
        """
        This generates slug for a project.
        :return:Returns slug.
        """
        # print "----->3"

        slug = self.project_title.replace(" ","-").replace("_","-")\
                .replace("/","-").replace("\\","-").replace(".","-").replace(",","-").replace(";",'-')\
                .replace(":","-").replace("--","-")


        # if check db if exist:

        slug_exists = db_object.do_select("select id from projects where slug = ?", (slug, )).fetchone()

        if slug_exists is not None:
            # TODO - adding suffix "-1" should not be manual
            slug = "%s-1" % slug
            return slug
        else:
            return slug
