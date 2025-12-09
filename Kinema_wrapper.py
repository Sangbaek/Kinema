import subprocess
import argparse
import numpy as np

def run_kinema_proton(E_beam, Q2, W, thetapq):
  try:
    command = ["./kinema", "-P", "-e", "{:.4f}".format(E_beam), "-Q", "{:.4f}".format(Q2), "-W", "{:.4f}".format(W), "-C", "{:.4f}".format(thetapq)]
    stdout, stderr = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    stdout = stdout.decode("utf-8")
    stdout = stdout.splitlines()

    for ind, line in enumerate(stdout):
      if ind == 5:
        theta_e = float(line[21:28])
      if ind == 6:
        p_e     = float(line[21:28])
      if ind == 11:
        theta_q = float(line[12:18])
      if ind == 34:
        proton_line = line
    p_p           = float(proton_line.split()[3])
    theta_pq_lab  = float(proton_line.split()[4])
    theta_pq_CM   = float(proton_line.split()[-2])
    print("{:.2f}, {:.2f}, {:.2f}, {:.2f}, 0".format(p_e, theta_e, p_p, theta_q - theta_pq_lab))  
    print("{:.2f}, {:.2f}, {:.2f}, {:.2f}, 180".format(p_e, theta_e, p_p, theta_q + theta_pq_lab))  

  except Exception as e:
    print(f"Error constructing command: {e}")
    return None

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Run Kinema for proton kinematics.")
  parser.add_argument("-e" , "--E_beam" , type=float, help="Beam energy in MeV")
  parser.add_argument("-Q"     , "--Q2"     , type=float, help="Momentum transfer squared in (MeV/c)^2")
  parser.add_argument("-W"      , "--W"      , type=float, help="Invariant mass in MeV")
  parser.add_argument("-C", "--thetapq", type=float, help="Angle between proton and momentum transfer in degrees")

  args = parser.parse_args()

  run_kinema_proton(args.E_beam, args.Q2, args.W, args.thetapq)