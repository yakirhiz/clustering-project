from invoke import task
import sys


@task(help={'k': '(int) number of clusters', "n": "(int) number of observations", "Random": "flag, add --Random to set Random == True"})
def run(c, k="no k", n="no n", Random=True):

    # Compile the C-API modules
    c.run("python3.8.5 km_setup.py build_ext --inplace")
    c.run("python3.8.5 sc_setup.py build_ext --inplace")

    # import main module
    import main

    if (k == "no k" or n == "no n") and (Random==False): # wrong input
        print("When \"Random\" Flag is disabled, the program must recieve k and n as inputs.")
        exit()
    elif (Random == True):
        main.main(0,0,True)
    else: # if Random is false, then 'k' and 'n' are provided
        try:
            int_k = int(k)
            int_n = int(n)
            
            assert(int_k > 0)
            assert(int_n > 0)
            assert(int_n > int_k)
            assert(isinstance(Random, bool))
            
            main.main(int_k, int_n, False)
        except:
            if not ((int(k) > 0) or (int(n) > 0) or (int(n) > k)):
                print("One of the follwing is a problem: (k > 0) or (n > 0) or (n > k).")    
            else:
                print("Problem with argument or something went worng.")
            exit()

    # remove compiled code
    c.run("rm *mySpecClust*.so")
    c.run("rm *kmeans*.so")


#debug
@task
def print_hello(c):
    print("Hellooooo!!!!!")

### Spectral Clustering tasks
@task
def build_sc(c):
    """runs >>python3.8.5 sc_setup.py build_ext --inplace"""
    c.run("python3.8.5 sc_setup.py build_ext --inplace")

@task(aliases=['del_sc'])
def delete_sc(c):
    c.run("rm *mySpecClust*.so")


@task
def rebuild_sc(c):
    c.run("rm -r build")
    c.run("rm *mySpecClust*.so")
    c.run("sleep 0.1")
    c.run("python3.8.5 sc_setup.py build_ext --inplace")


### Kmeans++ tasks
@task
def build_km(c):
    """runs >>python3.8.5 sc_setup.py build_ext --inplace"""
    c.run("python3.8.5 km_setup.py build_ext --inplace")

@task(aliases=['del_km'])
def delete_km(c):
    c.run("rm *kmeans*.so")


@task
def rebuild_km(c):
    c.run("rm -r build")
    c.run("rm *kmeans*.so")
    c.run("sleep 0.1")
    c.run("python3.8.5 km_setup.py build_ext --inplace")

