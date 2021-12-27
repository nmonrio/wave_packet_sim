import numpy as np
import csv

############################
print(np.exp(2))
print(np.exp([1, 2, 3]))
print(np.pi)

#############################

def save_psi_squared_as_csv(list_psi_squared):
    with open('psi squared.csv', 'w', newline='') as file:
        writer = csv.writer(file)

#        for i in list_psi_squared:
        writer.writerow(list_psi_squared)        
        file.close()

def psi_squared_2_txt(psi_squared):
        with open('psi_squared.txt', 'w') as file:
            for item in psi_squared:
                file.write(str(item) + "\n")

        file.close()

if __name__ == "__main__":
    psi_squared_2_txt([1, 2, 3, 4, 5])



