import calvados as cal
from argparse import ArgumentParser

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--config',nargs='?', default='config.yaml', const='config.yaml', type=str)
    args = parser.parse_args()

    fconfig = args.config

    cal.sim.run(fconfig=fconfig)
