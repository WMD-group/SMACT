from optparse import OptionParser

parser = OptionParser()
parser.add_option(
    "-i",
    "--IP",
    action="store",
    type="float",
    dest="IP",
    default=5.0,
    help="The first material's IP, default 5.0",
)
parser.add_option(
    "-e",
    "--EA",
    action="store",
    type="float",
    dest="EA",
    default=4.0,
    help="The first material's EA, default 4.0",
)
parser.add_option(
    "-w",
    "--window",
    action="store",
    type="float",
    dest="window",
    default=1.0,
    help="The window around the IP/EA to allow +/- , eq w=1 gives +/- 0.5. Default 1.0",
)
parser.add_option(
    "-g",
    "--gap",
    action="store",
    type="float",
    dest="gap",
    default=4.0,
    help="The bandgap above which a layer is considered insulating and disregarded Default 4.0",
)

(options, args) = parser.parse_args()


f = open("CollatedData.txt")
lines = f.readlines()
f.close()

HTL = []
ETL = []
conducting_ETL = []
conducting_HTL = []
window = options.window

for line in lines:
    inp = line.split()
    if inp[0] != "Species":
        Eg = float(inp[1])
        EA = float(inp[2])
        IP = float(inp[3])
        if Eg > 2.0:
            if EA >= options.EA - window * 0.5 and EA <= options.EA + window * 0.5:
                ETL.append(inp[0])
            if Eg < options.gap:
                conducting_ETL.append(inp[0])
        if IP <= options.IP + window * 0.5 and IP >= options.IP - window * 0.5:
            if EA < 3.9:
                HTL.append(inp[0])
                if Eg < options.gap:
                    conducting_HTL.append(inp[0])

print("Number of potential electron contacting layers: ", len(conducting_ETL))
print("Number of potential hole contacting layers: ", len(conducting_HTL))

print("Conductive electron contacting layers: ")
print(len(conducting_ETL))
print(conducting_ETL)
print("Conductive hole contacting layers: ")
print(len(conducting_HTL))
print(conducting_HTL)
