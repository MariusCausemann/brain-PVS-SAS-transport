from pvs_network_netflow import *

def test():

    print("")
    success = run_single_bifurcation_test()

    print("")
    success = run_three_junction_verification()

    print("")
    success = test_murray_data()

    print("")
    success = run_murray_tree()
    
if __name__ == "__main__":

    # FIXME: Update to include varying radius estimates. WRONG WITHOUT
    # IT. See Gjerde, Sanchez, Rognes, JAP (2023):
    # https://arxiv.org/pdf/2310.02429.pdf Section II and onwards.
    
    __DEBUG_MODE=False
    test()


