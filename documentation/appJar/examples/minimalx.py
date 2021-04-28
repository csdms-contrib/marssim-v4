import sys
import time
sys.path.append("../../")
from appJar import gui
def changeTab():
    pass
def do_documentation():
    pass
def do_edit():
    pass
def do_run():
    pass
def do_exit():
    pass
with gui("MARSSIM GUI",bg="blue") as app:
    app.size = (800,500)
    #app.addMenuList("Tabs", ["EXIT", "DOCUMENTATION", "EDIT PARAMA", "RUN MARSSIM"], changeTab)
    #app.setTabbedFrameDisableAllTabs("Tabs", False)
    with app.tabbedFrame("Tabs"):
        app.setTabbedFrameTabExpand("Tabs",expand=True)
        with app.tab("Documentation", bg="slategrey", sticky="new"):
            pass
            #with app.labelFrame("Login Form"):
                #pass
        with app.tab("Edit Parameters", sticky="ew", expand="all", bg="orangered"):
            pass
        with app.tab("Run Marssim"):
            pass
        with app.tab("Exit"):
            pass
            #with app.labelFrame("dnd", hideTitle=True, sticky="news", inPadding=(20,20)):
                #pass
        #with app.tab("Calculator", inPadding = (5,5)):
            #pass
            #app.label("calculator", "", bg="grey", relief="sunken", anchor="e")
            #pass
       # with app.tab("Panes"):
            #pass
            #with app.panedFrame("a", sticky = "news"):
            #pass
        #with app.tab("Labels", sticky = "news"):
            #pass
app.go()
