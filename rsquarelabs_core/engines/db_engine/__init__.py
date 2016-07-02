import sqlite3, os, logging
from rsquarelabs_core.config import  RSQ_DB_LOG, RSQ_DB_PATH

"""
db_engine is a module built as database handler across the project.
db_engine uses sqlite as the database .


"""


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
# create a file handler
handler = logging.FileHandler(RSQ_DB_LOG)
handler.setLevel(logging.INFO)
# create a logging format
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(handler)




class DBEngine:

    def __init__(self, db_name=RSQ_DB_PATH):
        self.do_connect(db_name)
        """
        self.conn , self.cur are instalised in do_connect, but we are installising this here now,
        so we might need to change do_connect(), but all we are doing there is instalising
        """
        self.conn = sqlite3.connect(db_name)
        self.cur = self.conn.cursor()





    def do_connect(self, db_name):
        """
        This provides connection to the database.
        :param db_name: Path to the database file
        :return: Returns connection object to the database
        """
        self.conn = sqlite3.connect(db_name)
        sql_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'schema.sql')
        tables_structure = open(sql_file_path).read()

        try:
            self.conn.executescript(tables_structure)
        except Exception as e:
            logger.debug(e)
            pass

        self.cur = self.conn.cursor()
        return self.conn

    def do_select(self, cmd, params):
        """
        This provides data which is selected from database using given parameters.
        :param cmd: command (query) as a string for selecting.
        :param params: Constraint parameters for selecting data.
        :return: The selected data returned as a tuple.
        """
        data = self.cur.execute(cmd, params)
        return data

    def do_insert(self, cmd, params):
        """
        This provides cursor for inserting the data into database using parameters.
        :param cmd: command (query) as a string for insertion.
        :param params: Constraint parameters for insertion data.
        :return: Returns cursor object which is inserted the data into database.
        """
        try:
            self.cur.execute(cmd, params)
            self.conn.commit()
        except Exception as e:
            logger.error(e)
            logger.error(e.message)
        return self.cur

    def do_update(self, cmd, params):
        """
        This provides cursor for updating the data into database using parameters.
        :param cmd: command (query) as a string for updating.
        :param params: Constraint parameters for updating data.
        :return: Returns cursor object which is updated the data into database.
        """
        try:
            self.cur.execute(cmd, params)
            self.conn.commit()
        except Exception as e:
            logger.error(e)
            logger.error(e.message)
        return self.cur