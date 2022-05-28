from tkinter import *

root = Tk()
root.geometry("800x600")
root.title("Calligraphy")

pen_up = Button(root, text="PEN UP")
pen_up.place(x=0, y=570, width=800, height=30)

canvas = Canvas(root, width=400, height=400)
canvas.pack()

canvas.create_oval(60,60,210,210,outline="black", fill="black")

# e = Entry(root, width=35, borderwidth=3)
# e.grid(row=0, column=0, columnspan=3, padx=10, pady=10)




# def myClick():
#     myLabel = Label(root, text="Clicked!")
#     myLabel.pack()
#

# myButton.pack()

root.mainloop()