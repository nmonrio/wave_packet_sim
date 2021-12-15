import os
import tkinter
from tkinter import *
from tkinter.filedialog import asksaveasfile
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)

# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from src.methods import wave_packet

import numpy as np


class App(Tk):
    def __init__(self, wave_packet):
        super().__init__()
        self.iconbitmap('../resources/icon.ico')
        self.title('Numerical Solution of 1D Schrödinger equation for a Gaussian wave packet')
        self.geometry('1200x700')
        self.minsize(1200, 700)
        self.wave_packet = wave_packet

        # Run button
        self.Button(self, text = 'Run simulation.')
        self.Button

    def set_initial_conds(self, wave_packet):
        """
        Method for setting the attributes of the wave function as inputs for the function.
        """
        for key in vars(wave_packet).keys():
            key = StringVar()
            key_entry = Entry(root, textvariable=key)
            key_entry.pack()
            key.trace("w", lambda *args: print(key.get()))

    def save(self):  # Save to path function
        Files = [('All Files', '*.*'),
                 ('Python Files', '*.py'),
                 ('Text Document', '*.txt')]
        file = asksaveasfile(filetypes=Files, defaultextension=Files)


# Save data as txt buttons




psi2_2_txt = Button(root, text='Psi Squared.', command=lambda: save())
psi2_2_txt.pack(side=TOP, pady=20)

tc_2_txt = Button(root, text='Transmission coefficient.', command=lambda: save())
tc_2_txt.pack(side=TOP, pady=20)

integral_2_txt = Button(root, text='Integral.', command=lambda: save())
integral_2_txt.pack(side=TOP, pady=20)

# Display tables
fig = Figure(figsize=(5, 4), dpi=200)
fig.add_subplot(111).plot(t, 2 * np.sin(2 * np.pi * t))
fig.add_subplot().plot()
fig.add_subplot().plot()
# Ideal plot here

canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
canvas.draw()
canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)

toolbar = NavigationToolbar2Tk(canvas, root)
toolbar.update()
canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)


def on_key_press(event):
    print("you pressed {}".format(event.key))
    key_press_handler(event, canvas, toolbar)


canvas.mpl_connect("key_press_event", on_key_press)


def _quit():
    root.quit()  # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
    # Fatal Python Error: PyEval_RestoreThread: NULL tstate


# Exit button
button = tkinter.Button(master=root, text="Click to exit.", command=_quit)
button.pack(side=tkinter.BOTTOM)

root.mainloop()

if __name__ == "__main__":
    root = App()
    root.mainloop()
