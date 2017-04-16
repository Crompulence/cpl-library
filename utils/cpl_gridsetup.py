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


from draw_grid import draw_grid

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

    def __init__(self, parent, initialvalues={}, **kwargs):
        wx.Panel.__init__(self, parent, **kwargs)

        self.space = wx.BoxSizer(wx.HORIZONTAL)
        self.space.Add(wx.StaticText(self,label="",size=(100,20)), 0)

        #Add number of cells input 
        self.nxyz = wx.Panel(self)
        self.nxyz.label = wx.StaticText(self,label="Cells",size=(50,-1))
        nx = initialvalues.get('nx',8)
        self.nxyz.nx = wx.SpinCtrl(self, value=str(nx))
        ny = initialvalues.get('ny',8)
        self.nxyz.ny = wx.SpinCtrl(self, value=str(ny))
        nz = initialvalues.get('nz',8)
        self.nxyz.nz = wx.SpinCtrl(self, value=str(nz))
        self.nx=self.nxyz.nx; self.ny=self.nxyz.ny; self.nz=self.nxyz.nz
        #Add to sizer
        self.nhbox = wx.BoxSizer(wx.HORIZONTAL)
        self.nhbox.Add(self.nxyz.label, 0)
        self.nhbox.Add(self.nxyz.nx, 0)
        self.nhbox.Add(self.nxyz.ny, 0)
        self.nhbox.Add(self.nxyz.nz, 0)
        #self.nxyz.SetSizer(self.nhbox)

        #Add number of processor input
        self.pxyz = wx.Panel(self)
        self.pxyz.label = wx.StaticText(self,label="Procs",size=(50,-1))
        px = initialvalues.get('px',1)
        self.pxyz.px = wx.SpinCtrl(self, value=str(px))
        py = initialvalues.get('py',1)
        self.pxyz.py = wx.SpinCtrl(self, value=str(py))
        pz = initialvalues.get('pz',1)
        self.pxyz.pz = wx.SpinCtrl(self, value=str(pz))
        self.px=self.pxyz.px; self.py=self.pxyz.py; self.pz=self.pxyz.pz
        #Add to sizer
        self.phbox = wx.BoxSizer(wx.HORIZONTAL)
        self.phbox.Add(self.pxyz.label, 0)
        self.phbox.Add(self.pxyz.px, 0)
        self.phbox.Add(self.pxyz.py, 0)
        self.phbox.Add(self.pxyz.pz, 0)
        #self.pxyz.SetSizer(self.phbox)

        #Add minimum domain size input
        self.min_xyz = wx.Panel(self)
        self.min_xyz.label = wx.StaticText(self,label="Min",size=(50,-1))
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
        #self.min_xyz.SetSizer(self.minhbox)

        #Add maximum domain size input
        self.max_xyz = wx.Panel(self)
        self.max_xyz.label = wx.StaticText(self,label="Max",size=(50,-1))
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
        #self.max_xyz.SetSizer(self.maxhbox)

        #Split panels with a vertical sizer
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.space, 0)
        self.vbox.Add(self.nhbox, 0)
        self.vbox.Add(self.phbox, 0)
        self.vbox.Add(self.minhbox, 0)
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

        #Add to sizer
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.Add(self.btn, 0)
        self.hbox.Add(self.save_b, 0)
        self.hbox.Add(self.gs, 0)
        self.SetSizer(self.hbox, 0)

class CntrlPanel(wx.Panel):

    def __init__(self, parent, **kwargs):
        wx.Panel.__init__(self, parent, **kwargs)

        #Add three panels
        self.CFDpanel = GridPanel(self,initialvalues={"ymin":0.75,"ymax":1.75})
        self.MDpanel = GridPanel(self)
        self.updatepanel = UpdatePanel(self)

        #Split panels with a vertical size
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.CFDpanel, 0)
        self.vbox.Add(self.MDpanel, 0)
        self.vbox.Add(self.updatepanel, 1)

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

    def update_grid(self, e):

        #Reset the figure ready for plot
        self.pyplotp.reset_fig()
        
        #CFD
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

        print("CFD",nx,ny,nz,xmin,ymin,zmin,xmax,ymax,zmax)          
        draw_grid(self.pyplotp.ax,nx,ny,nz,
                  px=px,py=py,pz=pz,
                  xmin=xmin,ymin=ymin,zmin=zmin,
                  xmax=xmax,ymax=ymax,zmax=zmax,
                  lc='r',fc='k',label="CFD",
                  XKCD_plots=self.XKCD_plots)

        #MD
        nx = self.cntrp.MDpanel.nx.GetValue()
        ny = self.cntrp.MDpanel.ny.GetValue()
        nz = self.cntrp.MDpanel.nz.GetValue()
        px = self.cntrp.MDpanel.px.GetValue()
        py = self.cntrp.MDpanel.py.GetValue()
        pz = self.cntrp.MDpanel.pz.GetValue()
        xmin = self.cntrp.MDpanel.xmin.GetValue()
        ymin = self.cntrp.MDpanel.ymin.GetValue()
        zmin = self.cntrp.MDpanel.zmin.GetValue()
        xmax = self.cntrp.MDpanel.xmax.GetValue()
        ymax = self.cntrp.MDpanel.ymax.GetValue()
        zmax = self.cntrp.MDpanel.zmax.GetValue()

        print("MD",nx,ny,nz,xmin,ymin,zmin,xmax,ymax,zmax)          
        draw_grid(self.pyplotp.ax,nx,ny,nz,
                  px=px,py=py,pz=pz,
                  xmin=xmin,ymin=ymin,zmin=zmin,
                  xmax=xmax,ymax=ymax,zmax=zmax,
                  fc='b',lc='g',label="MD",
                  XKCD_plots=self.XKCD_plots)

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
                print('Saving figure as ' + fpath)
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
             
def main():
    
    ex = wx.App()
    frame=Coupled_Grid(None)
    frame.Show(True)
    ex.MainLoop()    

if __name__ == '__main__':
    main()          
