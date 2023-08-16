# script contributed by Philippe Garteiser; garteiserp@omrf.ouhsc.edu
from pymol import cmd

aa_dicts = {
    "Clustal": {
        "hidrophobic": ("blue", "ala,ile,leu,met,phe,trp,val"),
        "positive"   : ("red", "lys,arg"),
        "negative"   : ("magenta", "glu,asp"),
        "polar"      : ("green", "asn,gln,ser,thr"),
        "cysteines"  : ("pink", "cys"),
        "glycines"   : ("orange", "gly"),
        "prolines"   : ("yellow", "pro"),
        "aromatic"   : ("cyan", "his,tyr"),
    },
    "Zappo": {
        "hidrophobic"             : ("0xFFAFAF", "ile,leu,val,ala,met"),
        "aromatic"                : ("0xFFC800", "phe,trp,tyr"),
        "positive"                : ("0x6464FF", "lys,arg,his"),
        "negative"                : ("0xFF0000", "asp,glu"),
        "hydrophilic"             : ("0x00FF00", "ser,thr,asn,gln"),
        "conformationally_special": ("0xFF00FF", "pro,gly"),
        "cysteine"                : ("0xFFFF00", "cys"),
    },
    "TurnPropensity": {
        "asn_max": ("0xFF0000", "asn"),
        "gly"      : ("0xFF0000", "gly"),
        "pro"      : ("0xF60909", "pro"),
        "asx"      : ("0xF30C0C", "asx"),
        "asp"      : ("0xE81717", "asp"),
        "ser"      : ("0xE11E1E", "ser"),
        "cys"      : ("0xA85757", "cys"),
        "tyr"      : ("0x9D6262", "tyr"),
        "lys"      : ("0x7E8181", "lys"),
        "xaa"      : ("0x7C8383", "xaa"),
        "gln"      : ("0x778888", "gln"),
        "trp"      : ("0x738C8C", "trp"),
        "thr"      : ("0x738C8C", "thr"),
        "arg"      : ("0x708F8F", "arg"),
        "his"      : ("0x708F8F", "his"),
        "glx"      : ("0x5BA4A4", "glx"),
        "glu"      : ("0x3FC0C0", "glu"),
        "ala"      : ("0x2CD3D3", "ala"),
        "phe"      : ("0x1EE1E1", "phe"),
        "met"      : ("0x1EE1E1", "met"),
        "leu"      : ("0x1CE3E3", "leu"),
        "val"      : ("0x07F8F8", "val"),
        "ile_min": ("0x00FFFF", "ile"),
    },
}

other = (
    "white",
    "not ala,ile,leu,met,phe,trp,val,lys,arg,glu,asp,asn,gln,ser,thr,cys,gly,pro,his,tyr",
)


def resicolor(palette="Clustal", selection="all"):
    """USAGE: resicolor <selection>
    colors all or the given selection with arbitrary
    coloring scheme.
    """
    for aaclass, color_aas in aa_dicts[palette].items():
        color, aas = color_aas

        cmd.select(f"{palette}_{aaclass}", f'resn {"+".join(aas.split())}')
        cmd.do(f"color {color}, {palette}_{aaclass} & {selection}")
    cmd.group(palette, f"{palette}*")

    cmd.hide("everything", "resn HOH")


cmd.extend("resicolor2", resicolor)
