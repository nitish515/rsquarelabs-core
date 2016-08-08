# Websuite Documentation


## Websuite

WebServer is built using a micro web framework - bottle (www.bottlepy.org). Its very lightweight with routing, sessions and other basic web needs which are good enough for this project at this moment.



## To create new route or page

new route means new page in webclient server(http://localhost:9090 ).

```python
# In server.py

@app.route('/websuite/myown-file.html')
def myown_file():
    print "Im new route"
    
    PROGRAMMING THE CODE HERE
    
    
    # get the file content from HTML_DIR(rsquarelabs-core/websuite/static/html/) 
    content =  open(os.path.join(HTML_DIR, 'myown-file.html')).read()
    # now = now helps to print the timestamp when the page is executed by server. 
    return template(content, now=now)    


@app.route('/websuite/project/:project_id')
def projects_view(project_id):
    print "Im accessing %s project" %project_id
    
    PROGRAMMING THE CODE HERE
    
    
    # get the file content from HTML_DIR(rsquarelabs-core/websuite/static/html/) 
    content =  open(os.path.join(HTML_DIR, 'project-view.html')).read()
    # now = now helps to print the timestamp when the page is executed by server. 
    return template(content, now=now)    

```
