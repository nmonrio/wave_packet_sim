from matplotlib import pyplot as plt
import pandas as pd
from matplotlib import animation
from numpy.core.fromnumeric import size


class results:
    """
    Object that implements the output methods for the wave packet.
    """

    def psi_squared_2_txt(psi_squared):
        """
        Outputs the values of the psi_squared array as a .txt file.
        """
        with open('psi_squared.txt', 'w') as file:
            for item in psi_squared:
                file.write(str(item) + "\n")

        file.close()

    def integral_2_txt(integral_values):
        """
        Outputs the values of the integral array as a .txt file.
        """
        with open('integral_results.txt', 'w') as file:
            for item in integral_values:
                file.write(str(item) + "\n")

        file.close()

    def trans_coeff_2_txt(trans_coef_values):
        """
        Outputs the values of the transmission coefficient array as a .txt file.
        """
        with open('trans_coeff.txt', 'w') as file:
            for item in trans_coef_values:
                file.write(str(item) + "\n")

        file.close()

    def psi_squared_2_plot(list_psi_squared, x, potential):
        """
        Method for plotting the probabilities.
        """
        plt.plot(x, list_psi_squared[0], color='black', label='Initial packet.')
        plt.plot(x, list_psi_squared[1], color='red', label='t = 500$\Delta$t')
        plt.plot(x, list_psi_squared[2], color='green', label='t = 1000$\Delta$t')
        plt.plot(x, list_psi_squared[3], color='blue', label='t = 1500$\Delta$t')
        plt.plot(x, list_psi_squared[4], color='orange', label='t = 2000$\Delta$t')
        plt.legend()
        plt.ylabel('$|\Psi(x)|^2$')
        plt.xlabel('$x$')
        plt.axis([-25, 10, 0, 0.5])
        vx = plt.twinx()
        vx.plot(x, potential, color='black')
        vx.set_ylim(0, 3)
        vx.set_ylabel('$V(x)$')

        # plt.savefig("psi_squared_plot.png")
        plt.show()

    def psi_squared_2_plot_ex3(list_psi_squared, x, potential):
        """
        Method for plotting the probabilities for exercise 3. Contains larger x and more timestamps
        """
        plt.plot(x, list_psi_squared[0], color='black', label='Initial packet.')
        plt.plot(x, list_psi_squared[1], color='red', label='t = 500$\Delta$t')
        plt.plot(x, list_psi_squared[2], color='green', label='t = 1500$\Delta$t')
        plt.plot(x, list_psi_squared[3], color='blue', label='t = 3000$\Delta$t')
        plt.plot(x, list_psi_squared[4], color='orange', label='t = 5000$\Delta$t')
        plt.plot(x, list_psi_squared[5], color='purple', label='t = 9000$\Delta$t')
        plt.legend()
        plt.ylabel('$|\Psi(x)|^2$')
        plt.xlabel('$x$')
        plt.axis([-100, 100, 0, 0.5])
        vx = plt.twinx()
        vx.plot(x, potential, color='black')
        vx.set_ylim(0, 3)
        vx.set_ylabel('$V(x)$')

        # plt.savefig("psi_squared_plot.png")
        plt.show()

    # def psi_squared_2_movie():

    def plot_trans_coefficient(trans_coeff_list, t_values, asymptote):
        plt.plot(t_values, trans_coeff_list, color='blue')
        plt.axhline(y=asymptote, color='black', linestyle='--')
        plt.ylabel('Transmission probability, $T(t)$ ')
        plt.xlabel('Time, t ($s$)')
        plt.axis([0, 60, 0, 0.004])
        plt.show()

    def plot_tc_log_df(filename):
        """
        Method for plotting the TC for different k0 using csv and dataframes.
        """

        # Since we are working with a lot of data, it is more efficient to use data frames and pandas library
        # instead of regular matplotlib.

        headers = ['Time, t (s)', '$k_0$ = 0.5', '$k_0$ = 1', '$k_0$ = 1.5', '$k_0$ = 2.0', '$k_0$ = 2.5',
                   '$k_0$ = 3.0',
                   '$k_0$ = 3.5', '$k_0$ = 4.0', '$k_0$ = 4.5', '$k_0$ = 5.0', '$k_0$ = 5.5', '$k_0$ = 6.0',
                   '$k_0$ = 6.5',
                   '$k_0$ = 7.0', '$k_0$ = 7.5', '$k_0$ = 8.0']
        df = pd.read_csv(filename, names=headers)
        df.set_index('Time, t (s)').plot()
        plt.ylabel('Transmission probability.')
        plt.yscale('log')
        plt.ylim(10 ** (-5), 1)

        plt.show()
