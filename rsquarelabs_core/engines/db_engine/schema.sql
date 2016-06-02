CREATE TABLE projects
       (id INTEGER PRIMARY KEY     AUTOINCREMENT,
       title          TEXT    NOT NULL,
       tags           TEXT    NOT NULL,
       user_email     TEXT    NOT NULL,
       short_note     TEXT    NOT NULL,
       slug           TEXT    NOT NULL,
       path           TEXT,
       config         TEXT,
       log            TEXT,
       type           TEXT,
       date           TEXT,
       is_delete      INT);
CREATE TABLE project_activity
    ( id  INTEGER PRIMARY KEY     AUTOINCREMENT,
    tool_name TEXT NOT NULL,
    step_no TEXT NOT NULL,
    step_name TEXT NOT NULL,
    command TEXT NOT NULL,
    pid TEXT,
    status TEXT NOT NULL,
    log_file TEXT NOT NULL,
    project_id INTEGER NOT NULL,
    created_at TEXT,
    updated_at TEXT);
CREATE TABLE project_files
    (id INTEGER PRIMARY KEY AUTOINCREMENT,
    file_name TEXT NOT NULL,
    file_content BLOB NOT NULL,
    project_id INTEGER NOT NULL,
    created_at TEXT,
    updated_at TEXT);
CREATE TABLE protocols
    (id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL,
    version INT,
    parent_protocol INT,
    master_protocol INT,
    protocol_data TEXT NOT NULL,
    is_delete INT)
