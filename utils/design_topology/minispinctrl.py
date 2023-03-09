'''
    MiniSpinCtrl.py

    A custom class that looks like the small gtk2 SpinCtrl, not the monstrous version under gtk3.
    Based on TextCtrl with integer checking
    An attempt has been made to make it work with python2 and python3 - See the variable PY2
    Because the idea behind this SpinCtrl is to allow it to be really small, size it to cope
    with your Min/Max values or use the SetFontSize function.
    Accepts manual input from number keys and the numeric keypad.

    Styles: wx.TE_RIGHT (Default)
            wx.TE_LEFT
            wx.TE_CENTRE

            wx.BORDER_NONE is always applied

    Events: EVT_MINISPIN        A value change occurred in the spinctrl
                                A value was keyed in
                                A mouse click up or down
                                A mouse scroll up or down
                                An arrow-up or arrow-down

            EVT_MINISPINUP      A click up
                                A mouse scroll up
                                An arrow-up key

            EVT_MINISPINDOWN    A click down
                                A mouse scroll down
                                An arrow-down key
    Event Notes:
        All events return GetValue() i.e.

            def OnSpin(self,event):
                current_spin_value = event.GetValue()

        If monitoring SpinUp/SpinDown and SpinCtrl, SpinUp or SpinDown is the first event, followed by
        the SpinCtrl

        Entering values manually gives a SpinCtrl event for each Valid key press
        i.e.
            typing in 12 will give an event for 1 and another for 2
            With 15 in the control, using delete or backspace to remove the 5 or the 1
            will only give an event after pressing Return
            With 11 in the control using arrow keys to position the cursor between the
            ones and then pressing 2 (121) will fire an event on pressing 2 but the returned
            value will be the number in the control 121

        Pasting a value in sends a single Spin event

        Entering values manually WHEN Limited is set to True
            If the value becomes less than the Minimum value, Minimum value is applied
            If the value becomes greater than the Maximum value, Maximum value is applied

    Functions:
        GetValue()              Returns numeric value in the control

        GetMin()                Returns minimum value

        GetMax()                Returns maximum value

        GetRange()              Returns tuple (min, max) values

        GetIncrement()          Get the value used to increment the control per spin

        GetFontSize()           Get the Font size used for the control

        SetValue(int)           Set numeric value in the control

        SetMin(int)             Set minimum value - Automatically resets Range

        SetMax(int)             Set maximum value - Automatically resets Range

        SetRange(min,max)       Set minimum and maximum values

        SetIncrement(int)       Set the value to increment the control per spin

        SetFontSize()           Set the Font size used for the control

        SetLimited(Bool)        Set whether the control should limit the value to fall
                                within the current range
                                If False values not within current range are displayed in red

        IsLimited()             Is the control currently limiting the value to fall
                                within the current range - Returns True/False

        SetBackgroundColour(colour)

        SetForegroundColour(colour)

        Enable(Bool)            Enable/Disable the control
                                On Disable the control value is frozen

        IsEnabled()             Is the control Enabled - Returns True/False

    Default Values:
        min     -       0
        max     -       100
        initial -       0
        style   -       Align Right | Border None
        limited -       True
        range   -       (min, max)
        increment       1
        font size       SYS_SYSTEM_FONT Size

Author:     J Healey
Created:    31/08/2018
Copyright:  J Healey - 2018
License:    GPL 2 or any later version
Email:      rolfofsaxony@gmail.com

Usage example:

import wx
import minispinctrl as MSC
class Frame(wx.Frame):

    def __init__(self, parent):
        wx.Frame.__init__ (self, parent, -1, "Mini Spin Control")

        panel = wx.Panel(self, -1, size=(200,200))
        #self.ctl = MSC.MiniSpinCtrl(panel, -1, pos=(-1,-1), size=(-1,-1), min=-5, max=15, initial=-5)
        self.ctl = MSC.MiniSpinCtrl(panel, -1, size=(50,20), min=-5, max=100, initial=0)
        self.ctl.SetBackgroundColour('gold')
        self.ctl.SetForegroundColour('green')
        self.ctl.Bind(MSC.EVT_MINISPIN, self.OnSpin)
        self.ctl.Bind(MSC.EVT_MINISPINUP, self.OnSpinUp)
        self.ctl.Bind(MSC.EVT_MINISPINDOWN, self.OnSpinDown)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.ctl, 0, 0, 0)
        self.SetSizer(sizer)
        self.Show()
        self.ctl.SetLimited(False)
        self.ctl.SetRange(-10,110)

    def OnSpin(self, event):
        obj = event.GetEventObject()
        print ("Spin", obj.GetValue())

    def OnSpinUp(self, event):
        print ("Spin Up", event.GetValue())

    def OnSpinDown(self, event):
        print ("Spin Down")

app = wx.App()
frame = Frame(None)
app.MainLoop()
'''

import wx
from wx.lib.embeddedimage import PyEmbeddedImage
import sys
if sys.version_info.major == 2:
    PY2 = True
else:
    PY2 = False
spinupdown = PyEmbeddedImage(
    b'iVBORw0KGgoAAAANSUhEUgAAADIAAABkCAYAAADE6GNbAAAABmJLR0QA/wD/AP+gvaeTAAAA'
    b'CXBIWXMAAAPiAAAD4gHuD5mHAAAAB3RJTUUH4ggSCB03FPpcVwAAAbZJREFUeNrtmrtKA0EY'
    b'Rs+gnZekSBQEQd9DMIUYL6go6qvYW4qILyUaExVBsLG1slIQxUvWJtNEjGMyszMr34GPrXaZ'
    b's7P7Mzv7gxBCCCGEEEJ8Y+8/SOwAGbBVZIndjoTNPGCKJGCAUpeETakoMgYo/yBhUy3KjGS/'
    b'5LwIElcOIm3gMmWJo84gM0eZwxQlNh0FurOR0su92KeEzWzsSuZSoVwTtSxPepKwGYsl0vIs'
    b'chJDoulZwqaRp8R+IAmbtTwk1gNL2KyGlFjJScKm7LuSGWACeMlZ5A2o+JSZ7lw0ixRv3ESU'
    b'yIBr18emFxPAVAJLoXvgQdsfQgghhBBCCCGEEEKIsLjsdM8B75HHeAs8Dnohl26GkPkAZnzd'
    b'lVYkiSc8/yOZAl4jiNSNMd6f0wrwnKPEUsiXrp6TxGkeVWQ5sMRZP4Ma6uOcu86xFuhGTedd'
    b'3xueZ6JNxOYan2V5hIgdQlVPEispLHHGgc8BJLZTWq8t9CnRTHEBukbCLU1/mxpjDv4gMpz6'
    b'p8GFQ5kdpSCt5b1EahSoP36MiN1xoWWOi/wpbTu0W/9hX2ALIYQQQgghhBDp8wVxiX7Ava8A'
    b'egAAAABJRU5ErkJggg==')
spindown = PyEmbeddedImage(
    b'iVBORw0KGgoAAAANSUhEUgAAADIAAABkCAYAAADE6GNbAAAABmJLR0QA/wD/AP+gvaeTAAAA'
    b'CXBIWXMAAAPiAAAD4gHuD5mHAAAAB3RJTUUH4ggSCBk1npj4fwAAAiJJREFUeNrt2M+LEmEc'
    b'BvDvK1tYEJg/Rr0EeulsRsHqevHgdPLmWbqFYtBJQpK5CN6yf8DSgwPhqUNERy8aaOcFL0ER'
    b'1cWMCEV5uqwgtKvt6rgjPB+Y28w77zPvO+983xEhIiIiIiIiIiIiIiIiIiIiIiIiIiIiIiKi'
    b'8wBgq76s649a04ASkaiI/D7j2qtWZxCRa0qp7roTD1Y9BaUU+v3+w2g0+ui0czRNE5fLZVmK'
    b'QCAg9Xr9iYh0T/qz2fQqlUpvT57OTg/TNF9ufZqbpvl5lyEajcZHS1768Xjs03X9yy5C6Lp+'
    b'DOC6ZYsOgGg8Hrc0xOHhIQDcsnw5nkwmD9xu99yKEJqm/QFwz/Llf9H4YDB4ZkWQXq+X2/lH'
    b'qtPpvHA4HFsJEAqF0O12yzsPMZvNBMBBLpd7t40g2Wz2NQDHpVUTAG5UKpXjiwZwOp0wDOMD'
    b'gI0rBMcm74tS6lcmk7mbTCZ/XKSNRCLxyev1JpRS00ut7RY3BxCLRCLnGo10Oj2bTqcR21XL'
    b'tVrtKBwO/1cIj8eD4XB4x3Yh6vW6iIi0Wq2nwWBwZQi/349ms/nYbluFf1Sr1VdnhfD5fBiN'
    b'Rs9tvxkrFApKRMQwjPfLAZRSEBGUy+U3S/sc++8qATjz+fzP5TDtdvsbgCu2n1IL8/l8EcYX'
    b'i8W+iwhSqdRXADf3JsQpI3NULBYB4P6iItjbnxcAbu/dSKwIw99SRERERERb8xeP3+hzsiuu'
    b'uQAAAABJRU5ErkJggg==')
spinup = PyEmbeddedImage(
    b'iVBORw0KGgoAAAANSUhEUgAAADIAAABkCAYAAADE6GNbAAAABmJLR0QA/wD/AP+gvaeTAAAA'
    b'CXBIWXMAAAPiAAAD4gHuD5mHAAAAB3RJTUUH4ggSCAYkOXLWEwAAAipJREFUeNrt2M+LEmEc'
    b'BvDvK9uyBQumqHgp9NLZhALdPQk6nfwT8q7YVUISj57KQK+WHhTDU4eI8KIXhTQ6Bgodimjx'
    b'YEKEMvJ0acCW3DZz/AHPBx68zAzz8H3fmUERIiIiIiJaAwC//e57iVt7W2ahxGkqlQKAuyIi'
    b'uq7vT4n5fG6UcASDwTMRQSQS+QLg+t5MZmESR4lE4puIwEij0fgK4MpelEkmk0pEJJvNvlks'
    b'oZSCiCCTybz8VUTt/FRyudzzxRKLcTgcGI/HT3b25kulkoiIVKvVh263G8uKiAhcLhcqlcqD'
    b'nV1i+Xz+1Ov1XljCiN1ux2AwuL2Lmzvo8/kuVcJINBrVZ7OZb2dKDAaD41AodPYvJYyEw+GP'
    b'hULhaB3LzLLqiUopAXBcr9ffNptNxyrXaLVaN0ejUQvAoVJbeJDpui4ADuLx+OtVJnE+sVjs'
    b'BQDLVjZ/u91+arFYsI4iHo8HnU4ns/F90e/3H62jwPl0u934xkpMp9N7NpttbkYRp9P5A8Ad'
    b'098xAPwnJycwo4SRQCAAADdMm8RkMnFomvbZzBJGNE37AOCaKZOp1WqfNlHCSLlcfrf2aaTT'
    b'6VebLGGkVqs9u+xUDi4qoZSSXq9X9Pv92p+OcTqdYrVaTduTxWLx/nA4fK+Uemzcz9IX9F8m'
    b'okTELyLfl5x7aPbDUkSuKqU6/72sduXjdK//iSEiIiIiIiIiIiIiIiIiIiIiIiIiIiIiItqe'
    b'nwRG71f6w806AAAAAElFTkSuQmCC')
spindisabled = PyEmbeddedImage(
    b'iVBORw0KGgoAAAANSUhEUgAAAJUAAAEqCAYAAAAcSRJbAAAABHNCSVQICAgIfAhkiAAAAAlw'
    b'SFlzAAAD4gAAA+IB7g+ZhwAAABl0RVh0U29mdHdhcmUAd3d3Lmlua3NjYXBlLm9yZ5vuPBoA'
    b'AADDSURBVHic7cExAQAAAMKg9U9tDQ+gAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
    b'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
    b'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
    b'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAPgytxAAARuD4uMAAAAASUVORK5C'
    b'YII=')

getspinupdownImage = spinupdown.GetImage
getspindownImage = spindown.GetImage
getspinupImage = spinup.GetImage
getspindisabledImage = spindisabled.GetImage

mscEVT_MINISPINCTRL = wx.NewEventType()
EVT_MINISPIN = wx.PyEventBinder(mscEVT_MINISPINCTRL, 1)
mscEVT_MINISPINUP = wx.NewEventType()
EVT_MINISPINUP = wx.PyEventBinder(mscEVT_MINISPINUP, 1)
mscEVT_MINISPINDOWN = wx.NewEventType()
EVT_MINISPINDOWN = wx.PyEventBinder(mscEVT_MINISPINDOWN, 1)

class SpinEvent(wx.PyCommandEvent):
    """ Event sent from the :class:`MiniSpinCtrl` when a spin value changes. """

    def __init__(self, eventType, eventId=1, value=0):
        """
        Default class constructor.

        :param `eventType`: the event type;
        :param `eventId`: the event identifier.
        """

        wx.PyCommandEvent.__init__(self, eventType, eventId)
        self._eventType = eventType
        self.value = value

    def GetValue(self):
        """
        Retrieve the value of the control at the time
        this event was generated."""
        return self.value

class MiniSpinCtrl(wx.Control):

    def __init__(self, parent, id=wx.ID_ANY, pos=wx.DefaultPosition, size=wx.DefaultSize, min=0, max=100, initial=0, style=wx.TE_RIGHT, name="MiniSpinCtrl"):
        """
        Default class constructor.

        @param parent:  Parent window. Must not be None.
        @param id:      identifier. A value of -1 indicates a default value.
        @param pos:     MiniSpinCtrl position. If the position (-1, -1) is specified
                        then a default position is chosen.
        @param size:    If the default size (-1, -1) is specified then a default size is chosen.
        @param min:     Minimum allowed value.
        @param max:     Maximum allowed value.
        @param initial: Initial value.
        @param style:   Alignment (Left,Middle,Right).
        @param name:    Widget name.
        """

        wx.Control.__init__(self, parent, id, pos=pos, size=size, name=name)
        self._min = min
        self._max = max
        self._initial = initial
        self._pos = pos
        self._size = size
        self._name = name
        self._limited = True
        self._id = id
        self._frozen_value = "0"
        self._fgcolour = wx.BLACK
        self._bgcolour = wx.WHITE
        self._outrange_colour = wx.RED
        self._increment = 1
        self._font = wx.SystemSettings.GetFont(wx.SYS_SYSTEM_FONT)

        if self._limited:
            if self._initial > self._max:
                self._initial = self._max
            if self._initial < self._min:
                self._initial = self._min

        self.SetForegroundColour(self._fgcolour)
        self.SetBackgroundColour(self._bgcolour)
        self._style = style | wx.BORDER_NONE
        # Initialize images
        self.InitialiseBitmaps()

        #MiniSpinCtrl

        self.ctl = wx.TextCtrl(self, id, value=str(self._initial), \
            pos=self._pos, size=self._size, style=self._style, name=self._name)
        if PY2:
            self.spinner = wx.StaticBitmap(self, -1, bitmap=wx.BitmapFromImage(self._img))
        else:
            self.spinner = wx.StaticBitmap(self, -1, bitmap=wx.Bitmap(self._img))

        #End

        # Bind the events
        self.Bind(wx.EVT_MOUSEWHEEL, self.OnScroll)
        self.spinner.Bind(wx.EVT_LEFT_DOWN, self.OnSpin)
        self.ctl.Bind(wx.EVT_CHAR_HOOK, self.OnChar)
        self.Bind(wx.EVT_PAINT, self.OnPaint)

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.ctl, 1, wx.EXPAND, 0)
        sizer.Add(self.spinner, 0, wx.EXPAND, 0)
        self.SetImage(self._initial)
        self.SetSizerAndFit(sizer)
        self.Show()

    def OnPaint(self,event):
        if self.IsInRange():
            self.ctl.SetForegroundColour(self.GetForegroundColour())
        else:
            self.ctl.SetForegroundColour(self._outrange_colour)
        self.ctl.SetBackgroundColour(self.GetBackgroundColour())
        self.spinner.SetBackgroundColour(self.GetBackgroundColour())

    def InitialiseBitmaps(self):
        self._imgup = self.SetImageSize(getspinupImage())
        self._imgdown = self.SetImageSize(getspindownImage())
        self._imgupdown = self.SetImageSize(getspinupdownImage())
        self._imgdisabled = self.SetImageSize(getspindisabledImage())
        if self._initial <= self._min:
            self._img = self._imgup
        elif self._initial >= self._max:
            self._img = self._imgdown
        else:
            self._img = self._imgupdown

    def SetImageSize(self,img):
        #Size the image as Full height and the width is half the height
        h = self.GetSize()[1]
        img = img.Scale(int(h/2),int(h),quality=wx.IMAGE_QUALITY_HIGH)
        return img

    def SetValue(self,value):
        self.ctl.SetValue(str(value))
        self.SetImage(value)
        self.Update()

    def GetValue(self):
        value = self.ctl.GetValue()
        try:
            value = int(value)
        except:
            value = 0
        return value

    def SetMin(self,value):
        self._min = value

    def GetMin(self):
        return self._min

    def SetMax(self,value):
        self._max = value

    def GetMax(self):
        return self._max

    def SetLimited(self,value):
        self._limited = value

    def SetRange(self,min,max):
        self._min = min
        self._max = max

    def GetRange(self):
        return self._min, self._max

    def IsLimited(self):
        return self._limited

    def IsEnabled(self):
        return wx.Control.IsEnabled(self)

    def Enable(self, value):
        if value and self.IsEnabled(): # If value = current state do nothing
            return
        if not value and not self.IsEnabled():
            return
        wx.Control.Enable(self, value)
        self.Update()
        if value:
            #Enable via callafter in case someone has been scrolling away on the disabled control
            wx.CallAfter(self.OnReset)
        else:
            #Disable - Freeze the controls Value and change bitmap
            self._frozen_value = self.ctl.GetValue()
            self._img = self._imgdisabled
            if PY2:
                self.spinner.SetBitmap(wx.BitmapFromImage(self._img))
            else:
                self.spinner.SetBitmap(wx.Bitmap(self._img))

    def OnReset(self):
        #Reset the control to the state it was in when it was Disabled
        self.ctl.SetValue(self._frozen_value)
        self.SetImage(int(self._frozen_value))

    def IsInRange(self):
        if self._min <= self.GetValue() <= self._max:
            return True
        else:
            return False

    def SetIncrement(self,value):
        self._increment = value

    def GetIncrement(self):
        return self._increment

    def SetFontSize(self,value):
        self._font.SetPointSize(value)
        self.ctl.SetFont(self._font)

    def GetFontSize(self):
        return self._font.GetPointSize()

    #Spin image clicked (Top half = Up | Bottom half = Down)
    def OnSpin(self, event):
        H = self._img.GetHeight() / 2
        pos = event.GetY()
        if pos < H:
            self.OnScroll(None, self._increment)
        else:
            self.OnScroll(None, -self._increment)

    #Keyboard input, test for Arrow Up / Arrow down
    def OnChar(self, event):
        obj = event.GetEventObject()
        pos = obj.GetInsertionPoint()
        curr_text = obj.GetValue()
        key = event.GetUnicodeKey()
        if key == wx.WXK_NONE:
            key = event.GetKeyCode()

        #Test for and Swap out numeric pad keys for standard 0-9 values
        if key == wx.WXK_NUMPAD0:
            key = 48
        elif key == wx.WXK_NUMPAD1:
            key = 49
        elif key == wx.WXK_NUMPAD2:
            key = 50
        elif key == wx.WXK_NUMPAD3:
            key = 51
        elif key == wx.WXK_NUMPAD4:
            key = 52
        elif key == wx.WXK_NUMPAD5:
            key = 53
        elif key == wx.WXK_NUMPAD6:
            key = 54
        elif key == wx.WXK_NUMPAD7:
            key = 55
        elif key == wx.WXK_NUMPAD8:
            key = 56
        elif key == wx.WXK_NUMPAD9:
            key = 57
        else:
            pass

        #Test for position keys
        if key == wx.WXK_UP:
            self.OnScroll(None, self._increment)
        elif key == wx.WXK_DOWN:
            self.OnScroll(None, -self._increment)
        elif key == wx.WXK_LEFT:
            event.Skip()
        elif key == wx.WXK_RIGHT:
            event.Skip()
        elif key == wx.WXK_BACK:
            event.Skip()
        elif key == wx.WXK_DELETE:
            event.Skip()
        elif key == wx.WXK_RETURN:
            new_text = curr_text
            #Test for beyond range and limited display
            if int(new_text) > self._max and self._limited:
                new_text = str(self._max)
            if int(new_text) < self._min and self._limited:
                new_text = str(self._min)
            self.ctl.SetValue(new_text)
            if len(new_text) > 0: #Avoid Return on empty text
                event = SpinEvent(mscEVT_MINISPINCTRL, self.GetId(),int(new_text))
                event.SetEventObject(self)
                self.GetEventHandler().ProcessEvent(event)
                self.SetImage(int(new_text))
            event.Skip()

        #Test for Minus character at first position
        elif chr(key) == "-" and pos == 0:
            event.Skip()

        #Test for Valid numeric input
        elif chr(key).isdigit():
            #Test for nothing in curr_text
            if len(curr_text) < 1:
                check_text = "0"
            else:
                check_text = curr_text

            #create replacement text comprising old text and the new character
            # Allowing for the insertion point being changed by the keys
            # Delete, Backspace and Left, Right arrow keys
            lsel = obj.GetSelection()[0]; rsel = obj.GetSelection()[1]
            selected = rsel - lsel
            if selected == 0:
                new_text = curr_text[:pos]+chr(key)+curr_text[pos:]
            #Remove whole string if selected
            elif selected == len(curr_text):
                new_text = chr(key)
            #Otherwise remove selected bit
            else:
                new_text = curr_text[:lsel]+chr(key)+curr_text[rsel:]

            #Test for beyond range and limited display
            if int(new_text) > self._max and self._limited:
                new_text = str(self._max)
            if int(new_text) < self._min and self._limited:
                new_text = str(self._min)

            self.ctl.SetValue(new_text)
            self.ctl.SetInsertionPoint(len(new_text))
            event = SpinEvent(mscEVT_MINISPINCTRL, self.GetId(),int(new_text))
            event.SetEventObject(self)
            self.GetEventHandler().ProcessEvent(event)
            self.SetImage(int(new_text))
            event.Skip()
        else:
            pass

    #Mouse scroll: check rotation for direction
    #Check for non event override value in spin
    def OnScroll(self, event=None, spin=None):
        value = self.ctl.GetValue()
        if value == "":
            value = "0"
        value = int(value)
        if event:
            if event.GetWheelRotation() > 0:
                rotation = self._increment
            else:
                rotation = -self._increment
        else:
            rotation = spin
        adj=True
        #All values are added, because negative rotation or spin are already negative values
        if rotation > 0:
            if value + rotation <= self._max and self._limited == True:
                value += rotation
            elif self._limited == False:
                value += rotation
            else:
                adj=False

            #fire SpinUp event if the adjustment value is not zero
            if adj:
                self.ctl.SetValue(str(value))
                event = SpinEvent(mscEVT_MINISPINUP, self.GetId(), value)
                event.SetEventObject(self)
                self.GetEventHandler().ProcessEvent(event)
                event = SpinEvent(mscEVT_MINISPINCTRL, self.GetId(), value)
                event.SetEventObject(self)
                self.GetEventHandler().ProcessEvent(event)

        elif rotation < 0:
            if value + rotation >= self._min and self._limited == True:
                value += rotation
            elif self._limited == False:
                value += rotation
            else:
                adj=False

            #fire SpinDown event if the adjustment value is not zero
            if adj:
                self.ctl.SetValue(str(value))
                event = SpinEvent(mscEVT_MINISPINDOWN, self.GetId(), value)
                event.SetEventObject(self)
                self.GetEventHandler().ProcessEvent(event)
                event = SpinEvent(mscEVT_MINISPINCTRL, self.GetId(), value)
                event.SetEventObject(self)
                self.GetEventHandler().ProcessEvent(event)
        else:
            #No rotation
            pass

        self.SetImage(value)

    def SetImage(self, value):
        #Set appropriate image
        if value <= self._min:
            self._img = self._imgup
        elif value >= self._max:
            self._img = self._imgdown
        else:
            self._img = self._imgupdown
        if PY2:
            self.spinner.SetBitmap(wx.BitmapFromImage(self._img))
        else:
            self.spinner.SetBitmap(wx.Bitmap(self._img))
        self.Layout()
