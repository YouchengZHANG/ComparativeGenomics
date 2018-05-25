

def plot(genomes):

    genomes = genomes
    interactome = {}

    fo = open(genomes + '.out','r+')
    f = fo.read().splitlines()

    for l in f:
        if l.startswith(genomes + '.'):
            if not l.split()[0] in interactome:
                interactome[l.split()[0]] = 1
            else:
                interactome[l.split()[0]] += 1

    size = len(interactome)
    degree = sum(interactome.values())
    connectivity = degree / size
    print('size =',size,'  total degree =',degree,'  average connectivity =',connectivity)

    #The x-axis of each plot presents the node degree for each node protein in the network and the other axis the frequency of each node degree.
    value = list(interactome.values())
    Z = {v : value.count(v) for v in value}
    #print(Z)

    plt.scatter(list(Z.keys()), list(Z.values()),alpha=0.5,marker='.')
    plt.show()

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    #plot('272561')
    plot('515635')