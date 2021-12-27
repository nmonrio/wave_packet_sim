import numpy as np
import csv
import matplotlib as plt

############################
#print(np.exp(2))
#print(np.exp([1, 2, 3]))
#print(np.pi)

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

def plot_tc_extended(tc_data_list, k_not_values, t_values):
    plt.plot(t_values, tc_data_list[0], label=f'$k_0$ ={k_not_values[0]}')
    plt.plot(t_values, tc_data_list[1], label=f'$k_0$ ={k_not_values[1]}')
    plt.plot(t_values, tc_data_list[2], label=f'$k_0$ ={k_not_values[2]}')
    plt.plot(t_values, tc_data_list[3], label=f'$k_0$ ={k_not_values[3]}')
    plt.plot(t_values, tc_data_list[4], label=f'$k_0$ ={k_not_values[4]}')
    plt.plot(t_values, tc_data_list[5], label=f'$k_0$ ={k_not_values[5]}')
    plt.plot(t_values, tc_data_list[6], label=f'$k_0$ ={k_not_values[6]}')
    plt.plot(t_values, tc_data_list[7], label=f'$k_0$ ={k_not_values[7]}')
    plt.plot(t_values, tc_data_list[8], label=f'$k_0$ ={k_not_values[8]}')
    plt.plot(t_values, tc_data_list[9], label=f'$k_0$ ={k_not_values[9]}')
    plt.plot(t_values, tc_data_list[10], label=f'$k_0$ ={k_not_values[10]}')
    plt.plot(t_values, tc_data_list[11], label=f'$k_0$ ={k_not_values[11]}')
    plt.plot(t_values, tc_data_list[12], label=f'$k_0$ ={k_not_values[12]}')
    plt.plot(t_values, tc_data_list[13], label=f'$k_0$ ={k_not_values[13]}')
    plt.plot(t_values, tc_data_list[14], label=f'$k_0$ ={k_not_values[14]}')
    plt.plot(t_values, tc_data_list[15], label=f'$k_0$ ={k_not_values[15]}')
    plt.ylabel('Transmission probability, $T(t)$ ')
    plt.xlabel('Time, t ($s$)')
    plt.axis([0, 60, 0, 0.004])
    plt.show()


if __name__ == "__main__":

    t_val = []
    for t in np.arange(0, 50, 0.005):
        t_val.append(t)

    psi_squared_2_txt(t_val)


