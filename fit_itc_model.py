#!/usr/bin/env python3

import argparse
import numpy
import scipy.optimize

# Command-line arguments
parser = argparse.ArgumentParser(
    description = 'Fit a binding model to ITC data.'
)
parser.add_argument(
    '-a', '--approximate_dilution',
    help = 'Use approximate expressions for dilution effects used by MicroCal',
    action = 'store_true',
    default = False
)
parser.add_argument(
    '-b', '--bootstrap_iterations',
    help = 'Number of iterations for bootstrapping confidence intervals.',
    type = int,
    default = 0
)
parser.add_argument(
    '-d', '--dissociation_width',
    help = 'Width of prior for log(dissociation constant) regularization term '
        'in units of k_B T. Default is one.',
    type = float,
    default = 1.0
)
parser.add_argument(
    '-e', '--enthalpy_width',
    help = 'Width of prior for binding enthalpy regularization term in units '
        'of k_B T. Default is one.',
    type = float,
    default = 1.0
)
parser.add_argument(
    '-g', '--guess',
    help = 'File containing initial guess for fit parameters.',
    type = str,
    default = None
)
parser.add_argument(
    '-i', '--independent_sites',
    help = 'Treat binding sites as independent and equivalent.',
    action = 'store_true',
    default = False
)
parser.add_argument(
    '-l', '--ligand_concentration',
    help = 'Concentration of ligand in syringe (micromolar). For multiple '
        'experiments, use a comma-separated list.',
    type = str,
    default = None
)
parser.add_argument(
    '-n', '--number_sites',
    help = 'Number of binding sites',
    type = int,
    default = 1
)
parser.add_argument(
    '-p', '--penalty',
    help = 'Penalty for regularization terms.',
    type = float,
    default = 0.0
)
parser.add_argument(
    '-r', '--receptor_concentration',
    help = 'Initial concentration of receptor in IIC cell (micromolar). For '
        'multiple experiments, use a comma-separated list.',
    type = str,
    default = None
)
parser.add_argument(
    '-s', '--skip',
    help = 'Number of injections to skip. For multiple experiments, use a '
        'comma-separated list.',
    type = str,
    default = '0'
)
parser.add_argument(
    '-t', '--temperature',
    help = 'Temperature in Kelvin. Used only for regularization prior widths.',
    type = float,
    default = 310.15
)
parser.add_argument(
    '-v', '--volume',
    help = 'Volume of ITC cell (microliters). For multiple experiments, use a '
        'comma-separated list.',
    type = str,
    default = None
)
parser.add_argument(
    '--print_cost',
    help = 'Print cost of starting guess without optimization.',
    action = 'store_true'
)
parser.add_argument(
    '--save_bootstrap',
    help = 'Write bootstrap samples of parameters and residuals to file.',
    type = str,
    default = None
)
parser.add_argument('itc_data_file',
    help = 'File(s) containing ITC data in the format from MicroCal PEAQ.',
    nargs = '+',
    type = str)
args = parser.parse_args()

# Calculate derived parameters and variances from fit parameters
# params = [DH_A1, DH_B1, DH_B2, KD_A1, KD_A2, KD_B2]
def get_twosite_DH_KD(params, variances):

    # Compute fourth binding enthalpy DH_A2 = DH_A1 + DH_B2 - DH_B1
    DH_A2 = params[0] + params[2] - params[1]
    DH_A2_var = variances[0] + variances[1] + variances[2]

    # Take logarithm of dissociation constants
    log_KD_A1 = numpy.log(params[3])
    log_KD_A2 = numpy.log(params[4])
    log_KD_B2 = numpy.log(params[5])
    log_KD_A1_var = variances[3] / numpy.square(params[3])
    log_KD_A2_var = variances[4] / numpy.square(params[4])
    log_KD_B2_var = variances[5] / numpy.square(params[5])

    # Compute log of fourth dissociation constant KD_B1 = KD_A1 * KD_B2 / KD_A2
    log_KD_B1 = log_KD_A1 + log_KD_B2 - log_KD_A2
    log_KD_B1_var = log_KD_A1_var + log_KD_B2_var + log_KD_A2_var

    # Compute microscopic cooperativity alpha = KD_A1 / KD_A2
    log_alpha = log_KD_A1 - log_KD_A2
    log_alpha_var = log_KD_A1_var + log_KD_A2_var

    # Compute second macroscopic dissociation constant
    # KD_2 = KD_A2 + KD_B2
    KD_2 = params[4] + params[5]
    log_KD_2 = numpy.log(KD_2)
    log_KD_2_var = (variances[4] + variances[5]) / numpy.square(KD_2)

    # Compute first macroscopic dissociation constant
    # KD_1 = 1 / (1 / KD_A1 + 1 / KD_B1) = KD_A1 * KD_B1 / (KD_A1 + KD_B1)
    # KD_1 = KD_A1 * KD_B2 / (KD_A2 + KD_B2) = KD_A1 * KD_B2 / KD_2
    log_KD_1 = log_KD_A1 + log_KD_B2 - log_KD_2
    log_KD_1_var = (log_KD_A1_var
        + (log_KD_A2_var + log_KD_B2_var) * numpy.square(params[4] / KD_2))

    # Compute macroscopic cooperativity gamma = 4 * KD_1 / KD_2
    # gamma = 4 * KD_A1 * KD_B2 / (KD_A2 + KD_B2)^2
    log_gamma = numpy.log(4) + log_KD_1 - log_KD_2
    log_gamma_var = (log_KD_A1_var
        + (numpy.square(params[4] - params[5]) * log_KD_B2_var
            + 4 * variances[4]) / numpy.square(KD_2))

    # Assert that site A is the site with higher affinity for first binding step
    if log_KD_A1 < log_KD_B1:

        DH = numpy.array([params[0], params[1], DH_A2, params[2]])
        DH_stdev = numpy.array([numpy.sqrt(v)
            for v in [variances[0], variances[1], DH_A2_var, variances[2]]])
        log_KD = numpy.array([log_KD_A1, log_KD_B1, log_KD_A2, log_KD_B2,
            log_alpha, log_KD_1, log_KD_2, log_gamma])
        log_KD_stdev = numpy.array([numpy.sqrt(v) for v in [
            log_KD_A1_var, log_KD_B1_var, log_KD_A2_var, log_KD_B2_var,
            log_alpha_var, log_KD_1_var, log_KD_2_var, log_gamma_var
        ]])

    else:

        DH = numpy.array([params[1], params[0], params[2], DH_A2])
        DH_stdev = numpy.array([numpy.sqrt(v) for v in [
            variances[1], variances[0], variances[2], DH_A2_var]])
        log_KD = numpy.array([log_KD_B1, log_KD_A1, log_KD_B2, log_KD_A2,
            log_alpha, log_KD_1, log_KD_2, log_gamma])
        log_KD_stdev = numpy.array([numpy.sqrt(v) for v in [
            log_KD_B1_var, log_KD_A1_var, log_KD_B2_var, log_KD_A2_var,
            log_alpha_var, log_KD_1_var, log_KD_2_var, log_gamma_var
        ]])

    return DH, DH_stdev, log_KD, log_KD_stdev

# Function to solve a cubic equation x^3 + a_2 * x^2 + a_1 * x + a_0 = 0
# Returns the root and its partial derivatives with respect to the coefficients
def cubic_solver(a2, a1, a0):

    # Define intermediate quantities
    X = a2 * a2 / 9 - a1 / 3

    # Decide on class of solution based on sign of X and magnitude of W
    if X < 0:

        sqrt_X = numpy.sqrt(-X)
        W_denom = 54 * X * sqrt_X
        W = ((9 * a1 - 2 * a2 * a2) * a2 - 27 * a0) / W_denom

        t = numpy.arcsinh(W) / 3
        trig_t = numpy.sinh(t)
        partial_W = sqrt_X * numpy.cosh(t) / numpy.sqrt(W * W + 1) / 1.5
        partial_X = -trig_t / sqrt_X

    else:

        sqrt_X = numpy.sqrt(X)
        W_denom = 54 * X * sqrt_X
        W = ((9 * a1 - 2 * a2 * a2) * a2 - 27 * a0) / W_denom

        if W > 1:

            t = numpy.arccosh(W) / 3
            trig_t = numpy.cosh(t)
            partial_W = sqrt_X * numpy.sinh(t) / numpy.sqrt(W * W - 1) / 1.5

        elif W < -1:

            t = numpy.arccosh(-W) / 3
            trig_t = -numpy.cosh(t)
            partial_W = sqrt_X * numpy.sinh(t) / numpy.sqrt(W * W - 1) / 1.5

        else:

            t = numpy.arccos(W) / 3
            trig_t = numpy.cos(t)

            if W == 1:
                partial_W = sqrt_X / 4.5
            else:
                partial_W = sqrt_X * numpy.sin(t) / numpy.sqrt(1 - W * W) / 1.5

        partial_X = trig_t / sqrt_X

    root = 2 * sqrt_X * trig_t - a2 / 3
    partial_a2 = partial_X * a2 / 4.5 + partial_W * (
        9 * a1 - 6 * a2 * a2 - 18 * a2 * W * sqrt_X) / W_denom - 1.0 / 3
    partial_a1 = partial_X / -3 + partial_W * (
        9 * a2 + 27 * W * sqrt_X) / W_denom
    partial_a0 = -27 * partial_W / W_denom

    return root, partial_a2, partial_a1, partial_a0

vec_cubic_solver = numpy.vectorize(cubic_solver)

# Class that concatenates residuals and Jacobian for multiple experiments
class MultiCost():

    def __init__(self, cost_list):

        self.cost_list = cost_list
        self.N_cost = len(cost_list)
        self.cost_indices = range(self.N_cost)
        self.jac_indices = [tuple(int(j >= i) for j in range(self.N_cost - 1))
            for i in self.cost_indices]

    def __call__(self, params, target_DH):

        return numpy.concatenate(
            [self.cost_list[i](
                params[numpy.r_[i, self.N_cost : params.size]], target_DH[i]
            ) for i in self.cost_indices]
        )

    def get_jac(self, params, target_DH):

        return numpy.concatenate(
            [numpy.insert(
                self.cost_list[i].get_jac(
                    params[numpy.r_[i, self.N_cost : params.size]], target_DH[i]
                ), self.jac_indices[i], 0, axis = 1
            ) for i in self.cost_indices]
        )

# Class that concatenates residuals and Jacobian for multiple experiments and
# then adds regularization terms
class MultiCostRegularized(MultiCost):

    def __init__(self, cost_list, reg_targets, reg_weights):

        # Call constructor of MultiCost base class
        super().__init__(cost_list)

        # Initialize parameters for regularization term
        self.reg_targets = reg_targets
        self.reg_weights = reg_weights

        # Jacobian of regularization term
        self.reg_jac = reg_weights * numpy.eye(reg_weights.size)

    def __call__(self, params, target_DH):

        residuals = super().__call__(params, target_DH)

        # Append regularization terms to residuals
        reg_params = numpy.concatenate((params[:-3], numpy.log(params[-3:])))
        return numpy.concatenate((residuals,
            self.reg_weights * (reg_params - self.reg_targets)))

    def get_jac(self, params, target_DH):

        residual_jac = super().get_jac(params, target_DH)

        # Regularization is applied to log(K_D), so J = 1 / K_D
        reg_param_jac = numpy.concatenate((numpy.ones(params.size - 3),
            params[-3:]))

        # Append Jacobian of regularization terms to Jacobian of residuals
        return numpy.concatenate((residual_jac, self.reg_jac / reg_param_jac))

# Binding model using approximate treatment of dilution as in MicroCal PEAQ
class IndependentSitesApprox():

    def __init__(self, skip, N_site, volume, RT, LT, curr_coeff, prev_coeff):

        # Initialize quantities not being fit
        self.skip = skip
        self.volume = volume
        self.NRT = N_site * RT
        self.LT = LT
        self.curr_coeff = curr_coeff
        self.prev_coeff = prev_coeff

        # Initialize Jacobian
        self.jacobian = numpy.zeros((LT.size, 4))

        # Jacobian of evolved heat with respect to offset
        self.jacobian[:, 0] = numpy.ones(self.jacobian.shape[0])

    def __call__(self, params, target_DH):

        # Parameters
        offset = params[0]
        eta = params[1]
        DH = params[2]
        KD = params[3]

        # Coefficient of linear term in quadratic expression for [L]
        b = eta * self.NRT - self.LT + KD

        # Square root of discriminant
        sqrt_discriminant = numpy.sqrt(b * b + 4 * self.LT * KD)

        # Concentration of free ligand [L] = sqrt(b^2 + 4 eta L_T K_D) - b
        free_ligand = (sqrt_discriminant - b) / 2

        # Jacobian of evolved heat with respect to binding enthalpy
        heat_jac_DH = self.volume * (self.LT - free_ligand)

        # Evolved heat Q = V_0 * DeltaH * (L_T - [L])
        heat = heat_jac_DH * DH

        # Jacobian of evolved heat with respect to eta
        heat_jac_eta = (self.volume * DH * free_ligand / sqrt_discriminant
            * self.NRT)

        # Jacobian of evolved heat with respect to dissociation constant
        heat_jac_KD = -heat / sqrt_discriminant

        # Residuals of injection enthalpy
        residuals = (self.curr_coeff * heat
            + numpy.insert(self.prev_coeff * heat[:-1], 0, 0) + offset
            - target_DH)

        # Jacobian of injetion enthalpy with respect to eta
        self.jacobian[:, 1] = ((self.curr_coeff * heat_jac_eta)
            + numpy.insert(self.prev_coeff * heat_jac_eta[:-1], 0, 0))

        # Jacobian of injetion enthalpy with respect to binding enthalpy
        self.jacobian[:, 2] = ((self.curr_coeff * heat_jac_DH)
            + numpy.insert(self.prev_coeff * heat_jac_DH[:-1], 0, 0))

        # Jacobian of injetion enthalpy with respect to dissociation constant
        self.jacobian[:, 3] = ((self.curr_coeff * heat_jac_KD)
            + numpy.insert(self.prev_coeff * heat_jac_KD[:-1], 0, 0))

        return residuals[self.skip:]

    def get_jac(self, params, target_DH):

        return self.jacobian[self.skip:]

# Binding model assuming independent and equivalent binding sites
class IndependentSites():

    def __init__(self, skip, N_site, R0, L0, RT, LT):

        # Initialize quantities not being fit
        self.skip = skip
        self.NRT = N_site * RT
        self.LT = LT
        self.NL0 = N_site * L0
        stoich = LT / RT
        self.stoich_diff = L0 * numpy.diff(stoich, prepend = 0)
        self.stoich_plus_stoich_0 = stoich + L0 / R0

        # Initialize Jacobian
        self.jacobian = numpy.zeros((LT.size, 4))

        # Jacobian of evolved heat with respect to offset
        self.jacobian[:, 0] = numpy.ones(self.jacobian.shape[0])

    def __call__(self, params, target_DH):

        # Parameters
        offset = params[0]
        eta = params[1]
        DH = params[2]
        KD = params[3]

        # Coefficient of linear term in quadratic expression for [L]
        b = eta * self.NRT - self.LT + KD

        # Square root of discriminant
        sqrt_discriminant = numpy.sqrt(b * b + 4 * self.LT * KD)

        # Concentration of free ligand [L] = sqrt(b^2 + 4 eta L_T K_D) - b
        free_ligand = (sqrt_discriminant - b) / 2

        # Integral of total enthalpy change with respect to stoichiometric ratio
        int_enthalpy = self.stoich_plus_stoich_0 * (self.LT - free_ligand)

        # Jacobian of injection enthalpy with respect to eta
        self.jacobian[:, 1] = (DH * self.NL0  / self.stoich_diff
            * numpy.diff(free_ligand / sqrt_discriminant, prepend = 0))

        # Jacobian of injection enthalpy with respect to binding enthalpy
        self.jacobian[:, 2] = (
            numpy.diff(int_enthalpy, prepend = 0) / self.stoich_diff)

        # Jacobian of injection enthalpy with respect to dissociation constant
        self.jacobian[:, 3] = (-DH / self.stoich_diff
            * numpy.diff(int_enthalpy / sqrt_discriminant, prepend = 0))

        # Residuals of injection enthalpy
        residuals = self.jacobian[:, 2] * DH + offset - target_DH

        return residuals[self.skip:]

    def get_jac(self, params, target_DH):

        return self.jacobian[self.skip:]

# Binding model for two interdependent, non-equivalent sites
class TwoSites():

    def __init__(self, skip, R0, L0, RT, LT):

        # Initialize quantities not being fit
        self.skip = skip
        self.RT = RT
        self.LT = LT
        self.L0 = L0
        stoich = LT / RT
        self.stoich_diff = L0 * numpy.diff(stoich, prepend = 0)
        self.stoich_plus_stoich_0 = stoich + L0 / R0

        # Initialize Jacobian
        self.jacobian = numpy.zeros((LT.size, 8))

        # Jacobian of evolved heat with respect to offset
        self.jacobian[:, 0] = numpy.ones(self.jacobian.shape[0])

    def __call__(self, params, target_DH):

        # Parameters
        offset = params[0]
        eta = params[1]
        DH_A1 = params[2]
        DH_B1 = params[3]
        DH_B2 = params[4]
        KD_A1 = params[5]
        KD_A2 = params[6]
        KD_B2 = params[7]

        eta_RT = eta * self.RT
        eta_RT_minus_LT = eta_RT - self.LT
        KD_A2_plus_KD_B2 = KD_A2 + KD_B2
        KD_A1_times_KD_B2 = KD_A1 * KD_B2

        # Coefficient of quadratic term in cubic expression for [L]
        # a2 = 2 * eta * RT - LT + KD_A2 + KD_B2
        a2 = eta_RT + eta_RT_minus_LT + KD_A2_plus_KD_B2

        # Coefficient of linear term in cubic expression for [L]
        # a1 = (eta * RT - LT) * (KD_A2 + KD_B2) + KD_A1 * KD_B2
        a1 = eta_RT_minus_LT * KD_A2_plus_KD_B2 + KD_A1_times_KD_B2

        # Coefficient of intercept term in cubic expression for [L]
        # a0 = -LT * KD_A1 * KD_B2
        a0 = -self.LT * KD_A1_times_KD_B2

        # Concentration of free ligand [L] and its partial derivatives with
        # respect to the cofficients a2, a1, and a0
        free_ligand, partial_free_a2, partial_free_a1, partial_free_a0 = (
            vec_cubic_solver(a2, a1, a0))

        # Concentration of bound ligand
        bound_ligand = self.LT - free_ligand

        # Normalization factor for concentrations of bound species
        bound_normalization = KD_A2_plus_KD_B2 + 2 * free_ligand

        # Common factor for terms in integratd enthalpy change
        int_common_factor = (self.stoich_plus_stoich_0 * bound_ligand
            / bound_normalization)

        # Bound species terms in integrated enthalpy change
        bound_species = (DH_A1 * KD_B2 + DH_B1 * KD_A2
            + (DH_A1 + DH_B2) * free_ligand)

        # Integral of total enthalpy change with respect to stoichiometric ratio
        int_enthalpy = int_common_factor * bound_species

        # Residuals of injection enthalpy
        residuals = (numpy.diff(int_enthalpy, prepend = 0) / self.stoich_diff
            + offset - target_DH)

        # Partial derivative of free ligand with respect to eta
        # partial_a2_eta = 2 * RT
        # partial_a1_eta = (KD_A2 + KD_B2) * RT
        # partial_a0_eta = 0
        partial_free_eta = self.RT * (
            2.0 * partial_free_a2 + KD_A2_plus_KD_B2 * partial_free_a1)

        # Partial derivative of free ligand with respect to dissociation
        # constants
        # partial_a2_KD_j = delta_{j,A2} + delta_{j,B2}
        # partial_a1_KD_j = (eta * RT - LT) * (delta_{j,A2} + delta_{j,B2})
        #     + KD_B2 * delta_{j,A1} + KD_A1 * delta_{j,B2}
        # partial_a0_KD_j = -LT (KD_B2 * delta_{j,A1} + KD_A1 * delta_{j,B2})
        partial_free_KD_A1 = KD_B2 * (
            partial_free_a1 - self.LT * partial_free_a0)
        partial_free_KD_A2 = (
            partial_free_a2 + eta_RT_minus_LT * partial_free_a1)
        partial_free_KD_B2 = (partial_free_a2
            + (eta_RT_minus_LT + KD_A1) * partial_free_a1
            - self.LT * KD_A1 * partial_free_a0)

        # Partial derivative of integrated enthalpy with respect to free ligand
        partial_int_DH_free_ligand = int_enthalpy * (
            (DH_A1 + DH_B2) / bound_species - 1.0 / bound_ligand
            - 2.0 / bound_normalization)

        # Jacobian of injection enthalpy with respect to eta
        self.jacobian[:, 1] = (
            numpy.diff(partial_int_DH_free_ligand * partial_free_eta,
                prepend = 0) / self.stoich_diff)

        # Jacobian of injection enthalpy with respect to binding enthalpies
        self.jacobian[:, 2] = (
            numpy.diff(int_common_factor * (KD_B2 + free_ligand),
                prepend = 0) / self.stoich_diff)
        self.jacobian[:, 3] = (
            numpy.diff(int_common_factor * KD_A2, prepend = 0)
            / self.stoich_diff)
        self.jacobian[:, 4] = (
            numpy.diff(int_common_factor * free_ligand, prepend = 0)
            / self.stoich_diff)

        # Jacobian of injection enthalpy with respect to dissociation constants
        normalized_bound_species = bound_species / bound_normalization
        self.jacobian[:, 5] = (
            numpy.diff(partial_int_DH_free_ligand * partial_free_KD_A1,
                prepend = 0) / self.stoich_diff)
        self.jacobian[:, 6] = (
            numpy.diff(partial_int_DH_free_ligand * partial_free_KD_A2
                + int_common_factor * (DH_B1 - normalized_bound_species),
                prepend = 0) / self.stoich_diff)
        self.jacobian[:, 7] = (
            numpy.diff(partial_int_DH_free_ligand * partial_free_KD_B2
                + int_common_factor * (DH_A1 - normalized_bound_species),
                prepend = 0) / self.stoich_diff)

        return residuals[self.skip:]

    def get_jac(self, params, target_DH):

        return self.jacobian[self.skip:]

# Binding model for two interdependent, non-equivalent sites with regularization
class TwoSitesRegularized(TwoSites):

    def __init__(self, skip, R0, L0, RT, LT, reg_targets, reg_weights):

        # Call constructor of TwoSites base class
        super().__init__(skip, R0, L0, RT, LT)

        # Initialize parameters for regularization term
        self.reg_targets = reg_targets
        self.reg_weights = reg_weights

        # Jacobian of regularization term
        self.reg_jac = reg_weights * numpy.eye(reg_weights.size)

    def __call__(self, params, target_DH):

        residuals = super().__call__(params, target_DH)

        # Append regularization terms to residuals
        reg_params = numpy.concatenate((params[:-3], numpy.log(params[-3:])))
        return numpy.concatenate((residuals,
            self.reg_weights * (reg_params - self.reg_targets)))

    def get_jac(self, params, target_DH):

        residual_jac = super().get_jac(params, target_DH)

        # Regularization is applied to log(K_D), so J = 1 / K_D
        reg_param_jac = numpy.concatenate((numpy.ones(params.size - 3),
            params[-3:]))

        # Append Jacobian of regularization terms to Jacobian of residuals
        return numpy.concatenate((residual_jac, self.reg_jac / reg_param_jac))

# Error checking on arguments
if args.ligand_concentration is None:

    print('You must provide the ligand concentration in micromolar (-l).')
    exit()

if args.receptor_concentration is None:

    print('You must provide the receptor concentration in micromolar (-r).')
    exit()

if args.volume is None:

    print('You must provide the cell volume in microliters (-v).')
    exit()

if args.penalty < 0:

    print('Regularization penalty (-p) must be non-negative.')
    exit()

if args.penalty > 0 and (args.approximate_dilution or args.independent_sites):

    print('Regularization penalty (-p) cannot be used with independent sites '
        '(-a or -i)')
    exit()

# Product of Boltzmann constant and temperature in kcal mol^-1
kBT = 0.0019872041 * args.temperature

# Get number of injections to skip, cell volume, and initial concentrations of
# ligand in syringe and receptor in cell
N_experiment = len(args.itc_data_file)
skip = numpy.array(
    numpy.ones(N_experiment) * numpy.array(args.skip.split(','), dtype = 'i8'),
    dtype = 'i8')
V0 = numpy.ones(N_experiment) * numpy.array(
    args.volume.split(','), dtype = 'f8')
L0 = numpy.ones(N_experiment) * numpy.array(
    args.ligand_concentration.split(','), dtype = 'f8')
R0 = numpy.ones(N_experiment) * numpy.array(
    args.receptor_concentration.split(','), dtype = 'f8')

# Load ITC data
delta_V = []
target_DH = []
RT = []
LT = []

for i in range(N_experiment):

    itc_data = numpy.loadtxt(args.itc_data_file[i], unpack = True)

    # Size of injections in microliters
    delta_V.append(itc_data[0])

    # Convert target enthalpy from cal mol^-1 to kcal mol^-1
    target_DH.append(itc_data[1] / 1000)

    # Total volume injected after each injection, i.e. cumulative sum of
    # delta_V, divided by the volume of the ITC cell
    fractional_volume = numpy.cumsum(itc_data[0]) / V0[i]

    # Get total concentration of receptor and ligand in cell after each
    # injection in micromolar
    if args.approximate_dilution:

        # Use the approximate expressions for dilution used in MicroCal PEAQ
        # R_T = R_0 * (1 - V / V_0 / 2) / (1 + V / V_0 / 2)
        # L_T = L_0 * (1 - V / V_0 / 2) * V / V_0
        RT.append(R0[i] * (2 - fractional_volume) / (2 + fractional_volume))
        LT.append(L0[i] * (1 - fractional_volume / 2) * fractional_volume)

    else:

        # Exact expressions for dilution
        # R_T = R_0 * exp(-V / V_0)
        # L_T = L_0 * (1 - exp(-V / V_0))
        exp_fractional_volume = numpy.exp(-fractional_volume)
        RT.append(R0[i] * exp_fractional_volume)
        LT.append(L0[i] * (1 - exp_fractional_volume))

# Decide on binding model
if args.approximate_dilution:

    # Approximate dilution. Requires independent and equivalent sites.
    # Need to fit offset, eta, enthalpy, and dissociation constant.

    # Coefficients of heat for current injection and for previous injection
    curr_coeff = [(0.5 / V0[i] + 1 / delta_V[i]) / L0[i]
        for i in range(N_experiment)]
    prev_coeff = [((0.5 / V0[i] - 1 / delta_V[i]) / L0[i])[1:]
        for i in range(N_experiment)]

    # Initialize cost function object
    if N_experiment == 1:

        cost = IndependentSitesApprox(
            skip[0], args.number_sites, V0[0], RT[0], LT[0], curr_coeff[0],
            prev_coeff[0]
        )
        cost_args = target_DH[0]

    else:

        cost = MultiCost(
            [IndependentSitesApprox(
                skip[i], args.number_sites, V0[i], RT[i], LT[i], curr_coeff[i],
                prev_coeff[i]
            ) for i in range(N_experiment)]
        )
        cost_args = target_DH

    guess = numpy.insert(numpy.zeros(N_experiment + 2), N_experiment, 1.0)
    lower_bounds = numpy.concatenate((
        numpy.full(N_experiment, -numpy.inf),
        numpy.array([0.0, -numpy.inf, 0.0])
    ))

elif args.independent_sites or args.number_sites == 1:

    # Independent and equivalent sites.
    # Need to fit offset, eta, enthalpy, and dissociation constant.

    # Initialize cost function object
    if N_experiment == 1:

        cost = IndependentSites(
            skip[0], args.number_sites, R0[0], L0[0], RT[0], LT[0]
        )
        cost_args = target_DH[0]

    else:

        cost = MultiCost(
            [IndependentSites(
                skip[i], args.number_sites, R0[i], L0[i], RT[i], LT[i]
            ) for i in range(N_experiment)]
        )
        cost_args = target_DH

    guess = numpy.insert(numpy.zeros(N_experiment + 2), N_experiment, 1.0)
    lower_bounds = numpy.concatenate((
        numpy.full(N_experiment, -numpy.inf),
        numpy.array([0.0, -numpy.inf, 0.0])
    ))

elif args.number_sites == 2:

    # Two interdependent, non-equivalent sites.
    # Need to fit offset, eta, 3 enthalpies, and 3 dissociation constants.

    if args.guess is None:

        # Fit with independent sites to obtain a starting guess for parameters

        if N_experiment == 1:

            ind_cost = IndependentSites(
                skip[0], args.number_sites, R0[0], L0[0], RT[0], LT[0]
            )
            cost_args = target_DH[0]

        else:

            ind_cost = MultiCost(
                [IndependentSites(
                    skip[i], args.number_sites, R0[i], L0[i], RT[i], LT[i]
                ) for i in range(N_experiment)]
            )
            cost_args = target_DH

        ind_guess = numpy.insert(
            numpy.zeros(N_experiment + 2), N_experiment, 1.0)
        ind_lower_bounds = numpy.concatenate((
            numpy.full(N_experiment, -numpy.inf),
            numpy.array([0.0, -numpy.inf, 0.0])
        ))
        ind_upper_bounds = numpy.full(ind_guess.size, numpy.inf)
        ind_guess = numpy.minimum(
            numpy.maximum(ind_guess, ind_lower_bounds), ind_upper_bounds)

        # Do least-squares fit with independent sites
        ind_fit = scipy.optimize.least_squares(
            ind_cost, ind_guess, jac = ind_cost.get_jac,
            bounds = (ind_lower_bounds, ind_upper_bounds), x_scale = 'jac',
            args = (cost_args,)
        )

        guess = numpy.insert(ind_fit.x, -1,
            [ind_fit.x[-2], ind_fit.x[-2], ind_fit.x[-1], ind_fit.x[-1]])

    # Set up cost with regularization
    if args.penalty > 0:

        # Targets are zero for offsets, one for eta, and parameters from fit
        # with independent sites for enthalpies and log(dissociation constants)
        regularization_targets = numpy.concatenate((numpy.zeros(N_experiment),
            numpy.array([1.0]), guess[-6:-3], numpy.log(guess[-3:])))

        # Widths (standard deviation of Gaussian prior) are 1 kcal mol^-1 for
        # offset, 0.05 for eta, k_B T for enthalpies, and k_B T for
        # log(dissociation constants)
        prior_widths = numpy.concatenate((
            numpy.ones(N_experiment), numpy.array([0.05]),
            numpy.full(3, args.enthalpy_width * kBT),
            numpy.full(3, args.dissociation_width)
        ))

        # Weight for regularization terms. The weight is penalty / (width)^2,
        # but take the sqrt() here because this will be applied to the
        # residuals, which are then squared and summed to form the cost function
        regularization_weights = numpy.sqrt(args.penalty) / prior_widths

        if N_experiment == 1:

            cost = TwoSitesRegularized(skip[0], R0[0], L0[0], RT[0], LT[0],
                regularization_targets, regularization_weights)
            cost_args = target_DH[0]

        else:

            cost = MultiCostRegularized(
                [TwoSites(skip[i], R0[i], L0[i], RT[i], LT[i])
                    for i in range(N_experiment)],
                regularization_targets, regularization_weights
            )
            cost_args = target_DH

    # Set up cost without regularization
    else:

        if N_experiment == 1:

            cost = TwoSites(skip[0], R0[0], L0[0], RT[0], LT[0])
            cost_args = target_DH[0]

        else:

            cost = MultiCost(
                [TwoSites(skip[i], R0[i], L0[i], RT[i], LT[i])
                    for i in range(N_experiment)]
            )
            cost_args = target_DH

    lower_bounds = numpy.concatenate((
        numpy.full(N_experiment, -numpy.inf),
        numpy.array([0.0, -numpy.inf, -numpy.inf, -numpy.inf, 0.0, 0.0, 0.0])
    ))

else:

    print('Models with more than two non-independent binding sites are not '
        'supported')
    exit()

# Initial guess for parameters
if args.guess == 'ones':
    guess = numpy.concatenate(
        (numpy.full(N_experiment, 0.0), numpy.full(7, 1.0)))
elif args.guess is not None:
    guess = numpy.loadtxt(args.guess)

# Boundaries on fit parameters
upper_bounds = numpy.full(guess.size, numpy.inf)
guess = numpy.minimum(numpy.maximum(guess, lower_bounds), upper_bounds)

# Cost of initial guess
if args.penalty > 0:
    guess_cost = numpy.sum(numpy.square(cost(guess, cost_args)[:-guess.size]))
else:
    guess_cost = numpy.sum(numpy.square(cost(guess, cost_args)))

if args.print_cost:

    print('%14.8f' % guess_cost)
    exit()

for i in range(N_experiment):

    print('# Experiment %d Ligand %14.8f uM Receptor %14.8f uM Volume %14.8f uL'
        % (i, L0[i], R0[i], V0[i]))

# Do least-squares fit
fit = scipy.optimize.least_squares(
    cost, guess, jac = cost.get_jac, bounds = (lower_bounds, upper_bounds),
    x_scale = 'jac', args = (cost_args,)
)

# Estimate goodness-of-fit by reduced chi square
# fit.cost is 0.5 * sum of square residuals, so multiply by 2
# If regularization is used, do not include regularization terms in cost
if args.penalty > 0:

    fit_cost = numpy.sum(numpy.square(fit.fun[:-fit.x.size]))
    chi_sq_dof = fit.fun.size - 2 * fit.x.size

else:

    fit_cost = 2 * fit.cost
    chi_sq_dof = fit.fun.size - fit.x.size

reduced_chi_sq = fit_cost / chi_sq_dof

# Estimate uncertainties in fit parameters
# Estimate covariance by Moore-Penrose pseudoinverse of (J^T J)
covariance = numpy.linalg.pinv(fit.jac.T.dot(fit.jac))
param_variance = numpy.diag(reduced_chi_sq * covariance)

# Estimate uncertainties in residuals by Gaussian propagation of uncertainty
fit_residual_stdev = numpy.sqrt(numpy.square(fit.jac).dot(param_variance))

# Confidence interval for residuals
residual_lower = numpy.array([
    fit.fun[i] - 1.96 * fit_residual_stdev[i] for i in range(fit.fun.size)])
residual_upper = numpy.array([
    fit.fun[i] + 1.96 * fit_residual_stdev[i] for i in range(fit.fun.size)])

# Get fit parameters and uncertainties
# Fit params is [offset, eta, enthalpy, free energy, dissociation constant]
if (args.approximate_dilution or args.independent_sites
    or args.number_sites == 1):

    log_KD = numpy.array([numpy.log(fit.x[-1])])
    log_KD_stdev = numpy.array([numpy.sqrt(param_variance[-1]) / fit.x[-1]])
    linear_params = numpy.concatenate((fit.x[:-1], kBT * log_KD))
    linear_stdev = numpy.concatenate((numpy.sqrt(param_variance[:-1]),
        kBT * log_KD_stdev))

    fit_params = numpy.append(linear_params, fit.x[-1])

else:

    DH, DH_stdev, log_KD, log_KD_stdev = get_twosite_DH_KD(
        fit.x[(N_experiment + 1):], param_variance[(N_experiment + 1):])
    linear_params = numpy.concatenate((fit.x[:(N_experiment + 1)], DH,
        kBT * log_KD[:4]))
    linear_stdev = numpy.concatenate((
        numpy.sqrt(param_variance[:(N_experiment + 1)]), DH_stdev,
        kBT * log_KD_stdev[:4]
    ))

    fit_params = numpy.concatenate((linear_params, numpy.exp(log_KD)))

# Construct confidence intervals [-1.96 stdev, +1.96 stdev]
confidence_lower = numpy.array(
    [linear_params[i] - 1.96 * linear_stdev[i]
        for i in range(len(linear_params))]
    + [numpy.exp(log_KD[i] - 1.96 * log_KD_stdev[i])
        for i in range(len(log_KD))]
)
confidence_upper = numpy.array(
    [linear_params[i] + 1.96 * linear_stdev[i]
        for i in range(len(linear_params))]
    + [numpy.exp(log_KD[i] + 1.96 * log_KD_stdev[i])
        for i in range(len(log_KD))]
)

print('# Guess_cost     %14.8f (kcal_mol^-1)^2 Guess_chi_sq     %14.8f' % (
    guess_cost, guess_cost / chi_sq_dof))
print('# Fit_cost       %14.8f (kcal_mol^-1)^2 Fit_chi_sq       %14.8f' % (
    fit_cost, reduced_chi_sq))

if args.bootstrap_iterations > 0:

    bootstrap_fit_params = numpy.zeros((args.bootstrap_iterations, guess.size))
    bootstrap_fit_residuals = numpy.zeros((args.bootstrap_iterations,
        fit.fun.size))

    for i in range(args.bootstrap_iterations):

        # Perturb target data by a Gaussian
        if N_experiment == 1:

            bootstrap_target_DH = numpy.concatenate((target_DH[0][:skip[0]],
                target_DH[0][skip[0]:] + fit.fun[:(target_DH[0].size - skip[0])]
                   * numpy.random.normal(0, 1, target_DH[0].size - skip[0])))

        else:

            bootstrap_target_DH = [numpy.concatenate((
                target_DH[i][:skip[i]],
                target_DH[i][skip[i]:] + fit.fun[:(target_DH[i].size - skip[i])]
                    * numpy.random.normal(0, 1, target_DH[i].size - skip[i])
            )) for i in range(N_experiment)]

        # Refit to perturbed target data
        bootstrap_fit = scipy.optimize.least_squares(
            cost, guess, jac = cost.get_jac,
            bounds = (lower_bounds, upper_bounds), x_scale = 'jac',
            args = (bootstrap_target_DH,)
        )

        bootstrap_fit_params[i] = bootstrap_fit.x
        bootstrap_fit_residuals[i] = bootstrap_fit.fun

    # Write bootstrap samples of residuals
    if args.save_bootstrap is not None:
        numpy.savetxt(args.save_bootstrap + '_residuals', bootstrap_fit_params)

    # Estimate fit parameters and confidence intervals
    bootstrap_median = numpy.median(bootstrap_fit_params, axis = 0)
    bootstrap_fit_lower = numpy.percentile(bootstrap_fit_params, 2.5, axis = 0)
    bootstrap_fit_upper = numpy.percentile(bootstrap_fit_params, 97.5, axis = 0)

    # Estimate cost from median of bootstrap residuals
    bootstrap_residuals = numpy.median(bootstrap_fit_residuals, axis = 0)
    if args.penalty > 0:

        bootstrap_cost = numpy.sum(numpy.square(
            bootstrap_residuals[:-fit.x.size]))

    else:

        bootstrap_cost = numpy.sum(numpy.square(bootstrap_residuals))

    # Estimate confidence interval for residuals from percentiles of bootstrap
    # fit residuals
    bootstrap_residual_lower = numpy.percentile(bootstrap_fit_residuals, 2.5,
        axis = 0)
    bootstrap_residual_upper = numpy.percentile(bootstrap_fit_residuals, 97.5,
        axis = 0)

    # Estimate bootstrap parameters and confidence intervals
    if (args.approximate_dilution or args.independent_sites
        or args.number_sites == 1):

        bootstrap_params = numpy.insert(bootstrap_median, -1,
            kBT * numpy.log(bootstrap_median[-1]))
        bootstrap_lower = numpy.insert(bootstrap_fit_lower, -1,
            kBT * numpy.log(bootstrap_fit_lower[-1]))
        bootstrap_upper = numpy.insert(bootstrap_fit_upper, -1,
            kBT * numpy.log(bootstrap_fit_upper[-1]))

        # Write bootstrap samples of parameters
        if args.save_bootstrap is not None:
            numpy.savetxt(args.save_bootstrap + '_params', bootstrap_fit_params)

    else:

        bootstrap_fit_var = numpy.var(bootstrap_fit_params, axis = 0, ddof = 1)
        derived_param_tuples = [
            get_twosite_DH_KD(bootstrap_fit_params[i, (N_experiment + 1):],
                bootstrap_fit_var[(N_experiment + 1):])
            for i in range(args.bootstrap_iterations)
        ]

        DH = numpy.array([t[0] for t in derived_param_tuples])
        log_KD = numpy.array([t[2] for t in derived_param_tuples])
        KD = numpy.exp(log_KD)

        bootstrap_params = numpy.concatenate((
            bootstrap_median[:(N_experiment + 1)], numpy.median(DH, axis = 0),
            kBT * numpy.median(log_KD[:, :4], axis = 0),
            numpy.median(KD, axis = 0)
        ))
        bootstrap_lower = numpy.concatenate((
            bootstrap_fit_lower[:(N_experiment + 1)],
            numpy.percentile(DH, 2.5, axis = 0),
            kBT * numpy.percentile(log_KD[:, :4], 2.5, axis = 0),
            numpy.percentile(KD, 2.5, axis = 0)
        ))
        bootstrap_upper = numpy.concatenate((
            bootstrap_fit_upper[:(N_experiment + 1)],
            numpy.percentile(DH, 97.5, axis = 0),
            kBT * numpy.percentile(log_KD[:, :4], 97.5, axis = 0),
            numpy.percentile(KD, 97.5, axis = 0)
        ))

        # Write bootstrap samples of parameters
        if args.save_bootstrap is not None:
            numpy.savetxt(args.save_bootstrap + '_params', numpy.concatenate((
                bootstrap_fit_params[:, :(N_experiment + 1)], DH, log_KD)))

    # Print bootstrap results
    def print_param(param_name, index):

        print(('# %-23s' + 6 * ' %14.8f') % (
            param_name, fit_params[index], confidence_lower[index],
            confidence_upper[index], bootstrap_params[index],
            bootstrap_lower[index], bootstrap_upper[index]
        ))

    print('# Bootstrap_cost %14.8f (kcal_mol^-1)^2 Bootstrap_chi_sq %14.8f' % (
        bootstrap_cost, bootstrap_cost / chi_sq_dof))
    print('# Parameter Fit_value Confidence_lower Confidence_upper '
        'Bootstrap_value Bootstrap_lower Bootstrap_upper')

    for i in range(N_experiment):

        print_param('Offset_%d_(kcal_mol^-1)' % (i + 1), i)

    print_param('Eta', N_experiment)

    if (args.approximate_dilution or args.independent_sites
        or args.number_sites == 1):

        print_param('DeltaH_(kcal_mol^-1)', N_experiment + 1)
        print_param('DeltaG_(kcal_mol^-1)', N_experiment + 2)
        print_param('KD_(uM)', N_experiment + 3)

    else:

        print_param('DeltaH_A1_(kcal_mol^-1)', N_experiment + 1)
        print_param('DeltaH_B1_(kcal_mol^-1)', N_experiment + 2)
        print_param('DeltaH_A2_(kcal_mol^-1)', N_experiment + 3)
        print_param('DeltaH_B2_(kcal_mol^-1)', N_experiment + 4)
        print_param('DeltaG_A1_(kcal_mol^-1)', N_experiment + 5)
        print_param('DeltaG_B1_(kcal_mol^-1)', N_experiment + 6)
        print_param('DeltaG_A2_(kcal_mol^-1)', N_experiment + 7)
        print_param('DeltaG_B2_(kcal_mol^-1)', N_experiment + 8)
        print_param('KD_A1_(uM)', N_experiment + 9)
        print_param('KD_B1_(uM)', N_experiment + 10)
        print_param('KD_A2_(uM)', N_experiment + 11)
        print_param('KD_B2_(uM)', N_experiment + 12)
        print_param('Alpha', N_experiment + 13)
        print_param('KD_1_(uM)', N_experiment + 14)
        print_param('KD_2_(uM)', N_experiment + 15)
        print_param('Gamma', N_experiment + 16)

    print('# Stoichiometric_ratio Target_enthalpy_(kcal_mol^-1) Fit_enthalpy '
        'Confidence_lower Confidence_upper Bootstrap_enthalpy Bootstrap_lower '
        'Bootstrap_upper')

    k = 0
    for i in range(len(target_DH)):

        k -= skip[i]
        for j in range(skip[i], target_DH[i].size):

            print('%14.8f %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f' % (
                LT[i][j] / RT[i][j], target_DH[i][j],
                fit.fun[k + j] + target_DH[i][j],
                residual_lower[k + j] + target_DH[i][j],
                residual_upper[k + j] + target_DH[i][j],
                bootstrap_residuals[k + j] + target_DH[i][j],
                bootstrap_residual_lower[k + j] + target_DH[i][j],
                bootstrap_residual_upper[k + j] + target_DH[i][j]
            ))

        k += target_DH[i].size

else:

    # Print fit results
    def print_param(param_name, index):

        print(('# %-23s' + 3 * ' %14.8f') % (
            param_name, fit_params[index], confidence_lower[index],
            confidence_upper[index]
        ))

    print('# Parameter Fit_value Confidence_lower Confidence_upper')

    for i in range(N_experiment):

        print_param('Offset_%d_(kcal_mol^-1)' % (i + 1), i)

    print_param('Eta', N_experiment)

    if (args.approximate_dilution or args.independent_sites
        or args.number_sites == 1):

        print_param('DeltaH_(kcal_mol^-1)', N_experiment + 1)
        print_param('DeltaG_(kcal_mol^-1)', N_experiment + 2)
        print_param('KD_(uM)', N_experiment + 3)

    else:

        print_param('DeltaH_A1_(kcal_mol^-1)', N_experiment + 1)
        print_param('DeltaH_B1_(kcal_mol^-1)', N_experiment + 2)
        print_param('DeltaH_A2_(kcal_mol^-1)', N_experiment + 3)
        print_param('DeltaH_B2_(kcal_mol^-1)', N_experiment + 4)
        print_param('DeltaG_A1_(kcal_mol^-1)', N_experiment + 5)
        print_param('DeltaG_B1_(kcal_mol^-1)', N_experiment + 6)
        print_param('DeltaG_A2_(kcal_mol^-1)', N_experiment + 7)
        print_param('DeltaG_B2_(kcal_mol^-1)', N_experiment + 8)
        print_param('KD_A1_(uM)', N_experiment + 9)
        print_param('KD_B1_(uM)', N_experiment + 10)
        print_param('KD_A2_(uM)', N_experiment + 11)
        print_param('KD_B2_(uM)', N_experiment + 12)
        print_param('Alpha', N_experiment + 13)
        print_param('KD_1_(uM)', N_experiment + 14)
        print_param('KD_2_(uM)', N_experiment + 15)
        print_param('Gamma', N_experiment + 16)

    print('# Stoichiometric_ratio Target_enthalpy_(kcal_mol^-1) Fit_enthalpy '
        'Confidence_lower Confidence_upper')

    k = 0
    for i in range(len(target_DH)):

        k -= skip[i]
        for j in range(skip[i], target_DH[i].size):

            print('%14.8f %14.8f %14.8f %14.8f %14.8f' % (
                LT[i][j] / RT[i][j], target_DH[i][j],
                fit.fun[k + j] + target_DH[i][j],
                residual_lower[k + j] + target_DH[i][j],
                residual_upper[k + j] + target_DH[i][j]
            ))

        k += target_DH[i].size

