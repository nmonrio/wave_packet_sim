from matplotlib import pyplot as plt
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
        plt.plot(x, list_psi_squared[1], color='red', label='t=500$\Delta$t')
        plt.plot(x, list_psi_squared[2], color='green', label='1000$\Delta$t')
        plt.plot(x, list_psi_squared[3], color='blue', label='1500$\Delta$t')
        plt.plot(x, list_psi_squared[4], color='orange', label='2000$\Delta$t')
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

    # def psi_squared_2_movie():

    def plot_trans_coefficient(trans_coeff_list, t_values, asymptote):
        plt.plot(t_values, trans_coeff_list, color='blue')
        plt.axhline(y=asymptote, color='black', linestyle='--')
        plt.ylabel('Transmission probability, $T(t)$ ')
        plt.xlabel('Time, t ($s$)')
        plt.axis([0, 60, 0, 0.004])
        plt.show()

    def log_plot_trans_coeff(trans_coeff_list, t_values, k_not):
        plt.plot(t_values, trans_coeff_list, color='blue', label='ko={}'.format(str(k_not)))
        plt.yscale('log')
        plt.legend()
        plt.ylabel('Transmission probability')
        plt.xlabel('Time')
        plt.axis([0, 50, 0.00001, 1.5])
        plt.show()
