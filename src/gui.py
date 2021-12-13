import os
import tkinter
from tkinter import *
from tkinter.filedialog import asksaveasfile

root = Tk()
root.iconbitmap('../resources/icon.ico')
root.title('Numerical Solution of 1D Schrödinger equation for a Gaussian wave packet')
root.geometry('1200x700')
root.minsize(1200, 700)

# Exercise holder
exercises_options = ['Exercise 1', 'Exercise 2', 'Exercise 3', 'Exercise 4', 'Exercise 5']
value_exercise = tkinter.StringVar(root)
value_exercise.set('Choose the exercise: ')
exercises_menu = tkinter.OptionMenu(root, value_exercise, *exercises_options)
exercises_menu.pack()


# Save data as txt buttons

def save():                                             # Save to path function
    Files = [('All Files', '*.*'),
             ('Python Files', '*.py'),
             ('Text Document', '*.txt')]
    file = asksaveasfile(filetypes=Files, defaultextension=Files)


psi2_2_txt = Button(root, text='Psi Squared.', command=lambda: save())
psi2_2_txt.pack(side=TOP, pady=20)

tc_2_txt = Button(root, text='Transmission coefficient.', command=lambda: save())
tc_2_txt.pack(side=TOP, pady=20)

integral_2_txt = Button(root, text='Integral.', command=lambda: save())
integral_2_txt.pack(side=TOP, pady=20)





# Code to add widgets will go here...
root.mainloop()
