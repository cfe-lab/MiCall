import Tkinter as tk
import ttk
import tkFileDialog

class MiCall:
    def __init__(self, master):
        self.__version__ = '0.1'
        self.master = master

root = tk.Tk()  # parent widget
root.wm_title('MiCall')

app = MiCall(root)

root.mainloop()  # enter Tk event loop
root.destroy()
