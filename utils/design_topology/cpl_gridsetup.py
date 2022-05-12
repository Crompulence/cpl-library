#!/usr/bin/python
# -*- coding: utf-8 -*-
import wx
from wx.lib.masked import NumCtrl
import wx.lib.agw.floatspin as FS

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.backends.backend_wxagg as wxaggb
import matplotlib.backends.backend_wx as wxb
import os
import signal
import subprocess as sp

from draw_grid import draw_grid
import inpututils
from latex2utf import latex2utf

class PyplotPanel(wx.Panel):

    def __init__(self, parent,**kwargs):
        wx.Panel.__init__(self,parent,**kwargs)
        self.parent = parent
        self.figure = matplotlib.figure.Figure()
        self.canvas = wxaggb.FigureCanvasWxAgg(self, -1, self.figure)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.Bind(wx.EVT_SIZE, self.sizeHandler)
        self.cmap = matplotlib.cm.RdYlBu_r

    def sizeHandler(self, event):
        self.canvas.SetSize(self.GetSize())

    def reset_fig(self):

        self.figure.clf(keep_observers=True)
        self.ax = self.figure.add_subplot(111)

    def redraw(self):
        #Set axis
        self.ax.set_xlim(self.ax.get_xlim())
        self.ax.set_ylim(self.ax.get_ylim())
        self.canvas.draw()

    def savefigure(self, fpath):
        fs = matplotlib.rcParams.get('font.size')
        matplotlib.rcParams.update({'font.size': 22})
        self.figure.savefig(str(fpath),dpi=300, transparent=True, 
                            bbox_inches='tight', pad_inches=0.1)
        matplotlib.rcParams.update({'font.size': fs})
        self.canvas.draw()


class GridPanel(wx.Panel):

    def __init__(self, parent, initialvalues={}, title="", cells=True, 
                 procs=True, Min=True, Max=True, minmaxcell=False, **kwargs):
        wx.Panel.__init__(self, parent, **kwargs)

        #minmax version of cells or proces
        assert (cells or procs) != minmaxcell

        #Add title
        self.title = wx.StaticText(self, label=title, 
                                   style= wx.ALIGN_CENTER | wx.SIMPLE_BORDER)
        font = wx.Font(18, wx.MODERN, wx.NORMAL, wx.NORMAL)
        self.title.SetFont(font)
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.Add(self.title)

        #Add number of cells input
        if cells or minmaxcell:
            self.nxyz = wx.Panel(self)
            if cells:
                self.nxyz.label = wx.StaticText(self,label="Cells",size=(50,-1))
                nx = initialvalues.get('nx',8)
                ny = initialvalues.get('ny',8)
                nz = initialvalues.get('nz',8)
            elif minmaxcell:
                self.nxyz.label = wx.StaticText(self,label="Mincell",size=(50,-1))
                nx = initialvalues.get('mincellx',1)
                ny = initialvalues.get('mincelly',1)
                nz = initialvalues.get('mincellz',1)

            self.nxyz.nx = wx.SpinCtrl(self, value=str(nx))
            self.nxyz.ny = wx.SpinCtrl(self, value=str(ny))
            self.nxyz.nz = wx.SpinCtrl(self, value=str(nz))
            self.nx=self.nxyz.nx; self.ny=self.nxyz.ny; self.nz=self.nxyz.nz
            #Add to sizer
            self.nhbox = wx.BoxSizer(wx.HORIZONTAL)
            self.nhbox.Add(self.nxyz.label, 0)
            self.nhbox.Add(self.nxyz.nx, 0)
            self.nhbox.Add(self.nxyz.ny, 0)
            self.nhbox.Add(self.nxyz.nz, 0)

        #Add number of processor input
        if procs or minmaxcell:
            self.pxyz = wx.Panel(self)
            if procs:
                self.pxyz.label = wx.StaticText(self,label="Procs",size=(50,-1))
                px = initialvalues.get('px',1)
                py = initialvalues.get('py',1)
                pz = initialvalues.get('pz',1)
            elif minmaxcell:
                self.pxyz.label = wx.StaticText(self,label="Maxcell",size=(50,-1))
                px = initialvalues.get('maxcellx',8)
                py = initialvalues.get('maxcelly',3)
                pz = initialvalues.get('maxcellz',8)

            self.pxyz.px = wx.SpinCtrl(self, value=str(px))
            self.pxyz.py = wx.SpinCtrl(self, value=str(py))
            self.pxyz.pz = wx.SpinCtrl(self, value=str(pz))
            self.px=self.pxyz.px; self.py=self.pxyz.py; self.pz=self.pxyz.pz
            #Add to sizer
            self.phbox = wx.BoxSizer(wx.HORIZONTAL)
            self.phbox.Add(self.pxyz.label, 0)
            self.phbox.Add(self.pxyz.px, 0)
            self.phbox.Add(self.pxyz.py, 0)
            self.phbox.Add(self.pxyz.pz, 0)

        #Add minimum domain size input
        if Min:
            self.min_xyz = wx.Panel(self)
            self.min_xyz.label = wx.StaticText(self,label="Origin",size=(50,-1))
            xmin = initialvalues.get('xmin',0.)
            self.min_xyz.xmin = FS.FloatSpin(self, -1, value=str(xmin))
            self.min_xyz.xmin.SetFormat("%f")
            self.min_xyz.xmin.SetDigits(2)
            ymin = initialvalues.get('ymin',0.)
            self.min_xyz.ymin = FS.FloatSpin(self, -1, value=str(ymin))
            self.min_xyz.ymin.SetFormat("%f")
            self.min_xyz.ymin.SetDigits(2)
            zmin = initialvalues.get('zmin',0.)
            self.min_xyz.zmin = FS.FloatSpin(self, -1, value=str(zmin))
            self.min_xyz.zmin.SetFormat("%f")
            self.min_xyz.zmin.SetDigits(2)
            self.xmin=self.min_xyz.xmin 
            self.ymin=self.min_xyz.ymin
            self.zmin=self.min_xyz.zmin
            #Add to sizer
            self.minhbox = wx.BoxSizer(wx.HORIZONTAL)
            self.minhbox.Add(self.min_xyz.label, 0)
            self.minhbox.Add(self.min_xyz.xmin, 0)
            self.minhbox.Add(self.min_xyz.ymin, 0)
            self.minhbox.Add(self.min_xyz.zmin, 0)

        #Add maximum domain size input
        if Max:
            self.max_xyz = wx.Panel(self)
            self.max_xyz.label = wx.StaticText(self,label="Domain",size=(50,-1))
            xmax = initialvalues.get('xmax',1.)
            self.max_xyz.xmax = FS.FloatSpin(self, -1, value=str(xmax))
            self.max_xyz.xmax.SetFormat("%f")
            self.max_xyz.xmax.SetDigits(2)
            ymax = initialvalues.get('ymax',1.)
            self.max_xyz.ymax = FS.FloatSpin(self, -1, value=str(ymax))
            self.max_xyz.ymax.SetFormat("%f")
            self.max_xyz.ymax.SetDigits(2)
            zmax = initialvalues.get('zmax',1.)
            self.max_xyz.zmax = FS.FloatSpin(self, -1, value=str(zmax))
            self.max_xyz.zmax.SetFormat("%f")
            self.max_xyz.zmax.SetDigits(2)
            self.xmax=self.max_xyz.xmax
            self.ymax=self.max_xyz.ymax
            self.zmax=self.max_xyz.zmax
            #Add to sizer
            self.maxhbox = wx.BoxSizer(wx.HORIZONTAL)
            self.maxhbox.Add(self.max_xyz.label, 0)
            self.maxhbox.Add(self.max_xyz.xmax, 0)
            self.maxhbox.Add(self.max_xyz.ymax, 0)
            self.maxhbox.Add(self.max_xyz.zmax, 0)

        #Split panels with a vertical sizer
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.hbox, 1, wx.ALL|wx.ALIGN_CENTER, 5)
        if cells or minmaxcell:
            self.vbox.Add(self.nhbox, 0)
        if procs or minmaxcell:
            self.vbox.Add(self.phbox, 0)
        if Min:
            self.vbox.Add(self.minhbox, 0)
        if Max:
            self.vbox.Add(self.maxhbox, 0)

        #When we set this sizer, causes a segfault
        self.SetAutoLayout(True)
        self.SetSizer(self.vbox)
        self.Layout()

class UpdatePanel(wx.Panel):

    def __init__(self, parent, initialvalues={}, **kwargs):
        wx.Panel.__init__(self, parent, **kwargs)

        #Bind button to regenerate grid
        self.btn = wx.Button(self, label='Update', pos=(-1, 250))

        #Add style button to grid
        self.gs = wx.CheckBox(self, label = 'Grid Style', pos = (-1,10))

        #Add savefig button
        self.save_b = wx.Button(self, label='SaveFig', pos=(-1, 250))

        self.run_b = wx.Button(self, label='RunCase', pos=(-1, 250))

        #Add to sizer
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.Add(self.btn, 0)
        self.hbox.Add(self.save_b, 0)
        self.hbox.Add(self.gs, 0)
        self.hbox.Add(self.run_b, 0)

        self.SetSizer(self.hbox, 0)

class TextPanel(wx.Panel):

    def __init__(self, parent, initialvalues={}, **kwargs):
        wx.Panel.__init__(self, parent, **kwargs)

        #Add some text
        font = wx.Font(18, wx.MODERN, wx.NORMAL, wx.NORMAL)
        import collections
        self.text = collections.OrderedDict()
        self.text["dx"] = wx.StaticText(self, wx.EXPAND, size=(200,-1))
        self.text["dy"] = wx.StaticText(self, wx.EXPAND, size=(200,-1))
        self.text["dz"] = wx.StaticText(self, wx.EXPAND, size=(200,-1))

        #Add to sizer
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        for t in self.text:
            self.text[t].SetFont(font)
            self.vbox.Add(self.text[t], 0)
        self.SetSizer(self.vbox, 0)

class CntrlPanel(wx.Panel):

    def __init__(self, parent, **kwargs):
        wx.Panel.__init__(self, parent, **kwargs)

        #Add three panels
        self.CFDpanel = GridPanel(self, initialvalues={"ymin":0.,"ymax":1.}, title="CFD")
        self.MDpanel = GridPanel(self, initialvalues={"nx":1,"ny":1,"nz":1}, title="MD", cells=False)
        self.CPLpanel = GridPanel(self, title="CPL", cells=False, procs=False, 
                                  Min=False, Max=False, minmaxcell=True)
        self.updatepanel = UpdatePanel(self)
        self.TextPanel = TextPanel(self)

        #Split panels with a vertical size
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.CFDpanel, 0)
        self.vbox.Add(self.MDpanel, 0)
        self.vbox.Add(self.CPLpanel, 0)
        self.vbox.Add(self.updatepanel, 1)
        self.vbox.Add(self.TextPanel, 1)

        self.SetAutoLayout(True)
        self.SetSizer(self.vbox)
        self.Layout()

class Coupled_Grid(wx.Frame):
           
    def __init__(self, *args, **kw):
        super(Coupled_Grid, self).__init__(*args, **kw) 
        
        self.XKCD_plots = False
        self.setup_panels()
        
    def setup_panels(self):

        #Instantiate plot panel and control panel
        self.pyplotp = PyplotPanel(self)
        self.cntrp = CntrlPanel(self)
    
        #Put into correct locations
        self.mainsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.mainsizer.Add(self.pyplotp, 1, wx.EXPAND | wx.ALL)
        self.mainsizer.Add(self.cntrp, 0, wx.EXPAND | wx.ALL)
        self.SetSizer(self.mainsizer)

        #Bind press of update button to plot update
        self.cntrp.updatepanel.btn.Bind(wx.EVT_BUTTON, self.update_grid)
        self.cntrp.updatepanel.gs.Bind(wx.EVT_CHECKBOX, self.onChecked)

        #Bind save figure button
        self.cntrp.updatepanel.save_b.Bind(wx.EVT_BUTTON, 
                                           lambda event: 
                                           self.save_dialogue(event, 'fig.png'))

        #Button to create and run case
        self.cntrp.updatepanel.run_b.Bind(wx.EVT_BUTTON, self.RunCase) 


    def update_grid(self, e):

        #Reset the figure ready for plot
        self.pyplotp.reset_fig()

        #CPL setup
        minnx = self.cntrp.CPLpanel.nx.GetValue()
        minny = self.cntrp.CPLpanel.ny.GetValue()
        minnz = self.cntrp.CPLpanel.nz.GetValue()
        maxnx = self.cntrp.CPLpanel.px.GetValue()
        maxny = self.cntrp.CPLpanel.py.GetValue()
        maxnz = self.cntrp.CPLpanel.pz.GetValue()

        print(("CPL", minnx, minny, minnz, maxnx, maxny, maxnz))          
        
        #MD
        px_MD = self.cntrp.MDpanel.px.GetValue()
        py_MD = self.cntrp.MDpanel.py.GetValue()
        pz_MD = self.cntrp.MDpanel.pz.GetValue()
        xmin_MD = self.cntrp.MDpanel.xmin.GetValue()
        ymin_MD = self.cntrp.MDpanel.ymin.GetValue()
        zmin_MD = self.cntrp.MDpanel.zmin.GetValue()
        xmax_MD = self.cntrp.MDpanel.xmax.GetValue()
        ymax_MD = self.cntrp.MDpanel.ymax.GetValue()
        zmax_MD = self.cntrp.MDpanel.zmax.GetValue()

        print(("MD",1,1,1,xmin_MD,ymin_MD,zmin_MD,xmax_MD,ymax_MD,zmax_MD))          
        draw_grid(self.pyplotp.ax,1,1,1,
                  px=px_MD,py=py_MD,pz=pz_MD,
                  xmin=xmin_MD,ymin=ymin_MD,zmin=zmin_MD,
                  xmax=xmax_MD,ymax=ymax_MD,zmax=zmax_MD,
                  fc='b',lc='g',label="MD",
                  XKCD_plots=self.XKCD_plots)

        #CFD
        nx_CFD = self.cntrp.CFDpanel.nx.GetValue()
        ny_CFD = self.cntrp.CFDpanel.ny.GetValue()
        nz_CFD = self.cntrp.CFDpanel.nz.GetValue()
        px_CFD = self.cntrp.CFDpanel.px.GetValue()
        py_CFD = self.cntrp.CFDpanel.py.GetValue()
        pz_CFD = self.cntrp.CFDpanel.pz.GetValue()
        xmin_CFD = self.cntrp.CFDpanel.xmin.GetValue()
        ymin_CFD = self.cntrp.CFDpanel.ymin.GetValue()
        zmin_CFD = self.cntrp.CFDpanel.zmin.GetValue()
        xmax_CFD = self.cntrp.CFDpanel.xmax.GetValue()
        ymax_CFD = self.cntrp.CFDpanel.ymax.GetValue()
        zmax_CFD = self.cntrp.CFDpanel.zmax.GetValue()
        dx = (xmax_CFD-xmin_CFD)/float(nx_CFD)
        dy = (ymax_CFD-ymin_CFD)/float(ny_CFD)
        dz = (zmax_CFD-zmin_CFD)/float(nz_CFD)

        ymin_CPL = ymax_MD - maxny*dy
        ymax_CPL = ymin_CPL + ny_CFD*dy

        print(("CFD",nx_CFD,ny_CFD,nz_CFD,xmin_CFD,ymin_CFD,zmin_CFD,xmax_CFD,ymax_CFD,zmax_CFD,dx,dy,dz))          
        draw_grid(self.pyplotp.ax,nx_CFD,ny_CFD,nz_CFD,
                  px=px_CFD,py=py_CFD,pz=pz_CFD,
                  xmin=xmin_CFD,ymin=ymin_CPL,zmin=zmin_CFD,
                  xmax=xmax_CFD,ymax=ymax_CPL,zmax=zmax_CFD,
                  lc='r',fc='k',label="CFD",
                  XKCD_plots=self.XKCD_plots)

        pnl = self.cntrp.TextPanel
        pnl.text["dx"].SetLabel(latex2utf("\Delta x=") + str(np.round(dx,5)))
        pnl.text["dy"].SetLabel(latex2utf("\Delta y=") + str(np.round(dy,5)))
        pnl.text["dz"].SetLabel(latex2utf("\Delta z=") + str(np.round(dz,5)))
        pnl.Layout()
        pnl.Update()
        #pnl.vbox.Layout()

        #Redraw
        self.pyplotp.redraw()

    def save_dialogue(self, event, defaultFile):
        dlg = wx.FileDialog(self, defaultDir='./', defaultFile=defaultFile,
                            style=wx.FD_SAVE) 
        if (dlg.ShowModal() == wx.ID_OK):
            fpath = dlg.GetPath()
        dlg.Destroy()

        #Check if defined, if cancel pressed then return
        try:
            fpath
        except NameError:
            return

        if defaultFile == 'fig.png':
            try:
                print(('Saving figure as ' + fpath))
                self.pyplotp.savefigure(fpath)
                print('Saved.')
            except ValueError:
                raise

    def onChecked(self, event):
        cb = event.GetEventObject()
        if cb.GetValue():
            self.XKCD_plots = True
        else:
            self.XKCD_plots = False

    def RunCase(self, e):

        #Reset the figure ready for plot
        self.pyplotp.reset_fig()

        #Get CFD inputs
        nx = self.cntrp.CFDpanel.nx.GetValue()
        ny = self.cntrp.CFDpanel.ny.GetValue()
        nz = self.cntrp.CFDpanel.nz.GetValue()
        px = self.cntrp.CFDpanel.px.GetValue()
        py = self.cntrp.CFDpanel.py.GetValue()
        pz = self.cntrp.CFDpanel.pz.GetValue()
        xmin = self.cntrp.CFDpanel.xmin.GetValue()
        ymin = self.cntrp.CFDpanel.ymin.GetValue()
        zmin = self.cntrp.CFDpanel.zmin.GetValue()
        xmax = self.cntrp.CFDpanel.xmax.GetValue()
        ymax = self.cntrp.CFDpanel.ymax.GetValue()
        zmax = self.cntrp.CFDpanel.zmax.GetValue()

        if any(a is 0 for a in [px,py,pz]):
            wx.MessageBox('Processor cannot be zero', 'Info', wx.OK)
            return
        
        #Change CFD input
        CFDDict = {"npxyz":[px, py, pz], 
                   "xyzL":[xmax, ymax, zmax], 
                   "xyz_orig":[xmin, ymin, zmin], 
                   "ncxyz":[nx, ny, nz]}
        ip_CFD = inpututils.InputMod("./CFD.in")
        for key in CFDDict:
            ip_CFD.replace_input(key, CFDDict[key])
        nproc_CFD = px*py*pz
        CFD_cmd = "mpiexec -n " + str(nproc_CFD) + " python ./CFD.py"

        px = self.cntrp.MDpanel.px.GetValue()
        py = self.cntrp.MDpanel.py.GetValue()
        pz = self.cntrp.MDpanel.pz.GetValue()
        xmin = self.cntrp.MDpanel.xmin.GetValue()
        ymin = self.cntrp.MDpanel.ymin.GetValue()
        zmin = self.cntrp.MDpanel.zmin.GetValue()
        xmax = self.cntrp.MDpanel.xmax.GetValue()
        ymax = self.cntrp.MDpanel.ymax.GetValue()
        zmax = self.cntrp.MDpanel.zmax.GetValue()

        if any(a is 0 for a in [px,py,pz]):
            wx.MessageBox('Processor cannot be zero', 'Info', wx.OK)
            return

        # Change MD input
        MDDict = {"npxyz":[px, py, pz], 
                  "xyzL":[xmax, ymax, zmax], 
                  "xyz_orig":[xmin, ymin, zmin]}
        ip_MD = inpututils.InputMod("./MD.in")
        for key in MDDict:
            ip_MD.replace_input(key, MDDict[key])
        nproc_MD = px*py*pz
        MD_cmd = "mpiexec -n " + str(nproc_MD) + " ./fortran/md"

        #CPL setup
        minnx = self.cntrp.CPLpanel.nx.GetValue()
        minny = self.cntrp.CPLpanel.ny.GetValue()
        minnz = self.cntrp.CPLpanel.nz.GetValue()
        maxnx = self.cntrp.CPLpanel.px.GetValue()
        maxny = self.cntrp.CPLpanel.py.GetValue()
        maxnz = self.cntrp.CPLpanel.pz.GetValue()

        #Change cpl input
        ip = inpututils.InputMod("./cpl/COUPLER.in")
        ip.replace_input("OVERLAP_EXTENTS", [minnx, maxnx, minny, maxny, minnz, maxnz])

        print(" =============== Running case =============== ")
        print("              ")
        print(MD_cmd)
        print("With")
        print(MDDict)
        print("              ")
        print(CFD_cmd)
        print("With")
        print(CFDDict)
        print("              ")
        print(" =============== Running case =============== ")
        cfd = sp.Popen("mpiexec -n " + str(nproc_CFD) + " python ./CFD.py", shell=True)
        fortran_md = "./fortran/md"
        if not os.path.isfile(fortran_md):
            os.chdir(fortran_md.replace("/md",""))
            sp.Popen("./build.sh", shell=True)
            os.chdir("../")
        print(fortran_md)
        md = sp.Popen("mpiexec -n " + str(nproc_MD) + " " + fortran_md , shell=True)

        #To ensure processes killed if this script is killed
        def kill_sp(pids):
            for pid in pids:
                if pid is None:
                    pass
                else:
                    try:
                        os.kill(pid, signal.SIGTERM)
                    except OSError:
                        pass
        import atexit
        atexit.register(kill_sp, [md.pid, cfd.pid])

        cfd.wait()
        md.wait()

#        from CFD import CFD
#        cfd = CFD()
#        cfd.recv_CPL_data()
#        cfd.plot_grid(self.pyplotp.ax)
#        cfd.plot_data(self.pyplotp.ax)
#        self.pyplotp.redraw()
#        cfd.finalise()
        
def main():
    
    ex = wx.App()
    frame=Coupled_Grid(None, size=(960,600))
    frame.Show(True)
    ex.MainLoop()    

if __name__ == '__main__':
    main()          
