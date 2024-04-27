from matplotlib import pyplot as plt


def main():
    map_rate = {
        'num': [0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0],
        '1A': [54.72, 53.64, 52.77, 52.03, 51.37, 50.77, 50.22, 49.76,
               49.22, 48.77, 48.34, 47.93, 47.53, 47.16, 46.80],
        '1B': [55.89, 54.67, 53.71, 52.90, 52.19, 51.55, 50.95, 50.47,
               49.91, 49.44, 48.99, 48.57, 48.16, 47.78, 47.41],
        '1D': [55.29, 53.98, 52.95, 52.06, 51.29, 50.60, 49.96, 49.43,
               48.83, 48.31, 47.83, 47.37, 46.93, 46.52, 46.12]
    }
    fig, ax = plt.subplots()
    ax.figure.set_size_inches(12, 8)
    ax.plot(map_rate['num'], map_rate['1A'], label='chr1A', alpha=0.8)
    ax.plot(map_rate['num'], map_rate['1B'], label='chr1B', alpha=0.8)
    ax.plot(map_rate['num'], map_rate['1D'], label='chr1D', alpha=0.8)
    ax.set_xlabel('Depth', fontsize=20)
    ax.set_ylabel('Mapping rate (%)', fontsize=20)
    ax.set_title('Mapping rate & Depth', fontsize=28)
    ax.legend(fontsize=18, framealpha=0.5)
    plt.show()


if __name__ == '__main__':
    main()
