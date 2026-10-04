import numpy as np

def clarke_transform(i_a: float, i_b: float) -> tuple:
    """
    Parameters:
    i_a : float : Phase A current
    i_b : float : Phase B current
    Returns:
    i_alpha : float : Alpha-axis current (directly aligned with Phase A)
    i_beta : float : Beta-axis current (perpendicular component)
    """
    # Reconstruct the third phase current
    i_c = -i_a - i_b  # From Kirchoff's Current Law

    # Clarke Transformation equations
    i_alpha = i_a  # Directly align with Phase A (amplitude-invariant)
    i_beta = (1/np.sqrt(3)) * (i_b - i_c)   # Orthogonal to alpha-axis


    return i_alpha, i_beta

def park_transform(i_alpha: float, i_beta: float, theta_e: float) -> tuple:
    """
    Parameters:
    i_alpha : float :Alpha-axis current (stationary frame)
    i_beta : float :Beta-axis current (stationary frame)
    theta_e : float :Rotor electrical angle in radians.
    Returns:
    i_d : float :d-axis current (rotor-aligned flux-related current)
    i_q : float :q-axis current (rotor-aligned torque-producing current)
    """
    # Park Transformation equations
    i_d =  i_alpha * np.cos(theta_e) + i_beta * np.sin(theta_e)
    i_q = -i_alpha * np.sin(theta_e) + i_beta * np.cos(theta_e)

    return i_d, i_q

def inverse_park(v_d: float, v_q: float, theta_e: float) -> tuple:
    """
    Parameters:
    v_d : float : Voltage aligned with the d-axis
    v_q : float : Voltage aligned with the q-axis
    theta_e : float : Rotor electrical angle in radians.
    Returns:
    v_alpha : float : Alpha-axis voltage (stationary frame)
    v_beta : float : Beta-axis voltage (stationary frame)
    """
    # Inverse Park Transformation equations
    v_alpha = v_d * np.cos(theta_e) - v_q * np.sin(theta_e)
    v_beta = v_d * np.sin(theta_e) + v_q * np.cos(theta_e)

    return v_alpha, v_beta

#test block 
"""if __name__ == "__main__":
    # Example inputs
    i_a = 10.0
    i_b = -5.0
    theta_e = np.pi / 4  # 45 degrees in radians

    # Clarke Transform Test
    i_alpha, i_beta = clarke_transform(i_a, i_b)
    print(f"Clarke Transform Test Results:")
    print(f"  i_alpha = {i_alpha:.3f}")
    print(f"  i_beta  = {i_beta:.3f}")

    # Park Transform Test
    i_d, i_q = park_transform(i_alpha, i_beta, theta_e)
    print(f"Park Transform Test Results:")
    print(f"  i_d = {i_d:.3f}")
    print(f"  i_q = {i_q:.3f}")

    # Inverse Park Transform Test
    v_d, v_q = 150.0, -50.0
    v_alpha, v_beta = inverse_park(v_d, v_q, theta_e)
    print(f"Inverse Park Transform Test Results:")
    print(f"  v_alpha = {v_alpha:.3f}")
    print(f"  v_beta  = {v_beta:.3f}")"""