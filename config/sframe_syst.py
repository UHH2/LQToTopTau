#!/usr/bin/env python

import os, sys, itertools

sys.path.append('/nfs/dust/cms/user/tholenhe/installs/varial-stable/Varial')


##################################### definition of UserConfig item changes ###

sys_uncerts = {
    # 'name' : {'item name': 'item value', ...},
    # 'TAUID_up'        : {'TauIdVariation':'up'},
    # 'TAUID_down'        : {'TauIdVariation':'down'},
    # 'TAUCH_up'        : {'TauChargeVariation':'up'},
    # 'TAUCH_down'        : {'TauChargeVariation':'down'},
    # 'SCALE_upup'        : {'ScaleVariationMuR':'up','ScaleVariationMuF':'up'},
    # 'SCALE_downdown'        : {'ScaleVariationMuR':'down','ScaleVariationMuF':'down'},
    # 'PU_up'        : {'puvariation':'up'},
    # 'PU_down'        : {'puvariation':'down'},
     #'MUID_up'        : {'muonvariation':'up'},
     #'MUID_down'        : {'muonvariation':'down'},
     #'MUTR_up'        : {'muontriggervariation':'up'},
     #'MUTR_down'        : {'muontriggervariation':'down'},
   #  'BTAG_up'        : {'OS_MCBtagEfficiencies':'up','SS_MCBtagEfficiencies':'up'},
   #  'BTAG_down'        : {'OS_MCBtagEfficiencies':'down','SS_MCBtagEfficiencies':'down'},
     #'MUISO_up'        : {'Systematic_MuonIso':'up'},
     #'MUISO_down'        : {'Systematic_MuonIso':'down'},
   #  'TAUEN_up'        : {'TauEnergyVariation':'up'},
   #  'TAUEN_down'        : {'TauEnergyVariation':'down'},
     'JEC_up'        : {'jecsmear_direction':'up'},
     'JEC_down'      : {'jecsmear_direction':'down'},
   #  'JER_up'        : {'jersmear_direction':'up'},
   #  'JER_down'      : {'jersmear_direction':'down'},
    # 'jer_jec_up'    : {'jersmear_direction':'up','jecsmear_direction':'up'},
    # 'jer_jec_down'  : {'jersmear_direction':'down','jecsmear_direction':'down'},
}
start_all_parallel = True

############################################################### script code ###
import varial
import sys
import os

if len(sys.argv) != 2:
    print 'Plz. give me da name of da sframe-config! ... dude!'
    exit(-1)


def set_uncert_func(uncert_name):
    uncert = sys_uncerts[uncert_name]
    def do_set_uncert(element_tree):
        cycle = element_tree.getroot().find('Cycle')
        user_config = cycle.find('UserConfig')
        output_dir = cycle.get('OutputDirectory')

        cycle.set('OutputDirectory', os.path.join(output_dir, uncert_name+'/'))

        for name, value in uncert.iteritems():
            uc_item = list(i for i in user_config if i.get('Name') == name)
            assert uc_item, 'could not find item with name: %s' % name
            uc_item[0].set('Value', value)

    return do_set_uncert


from varial.extensions.sframe import SFrame
from varial import tools
if start_all_parallel:
    ToolChain = tools.ToolChainParallel
else:
    ToolChain = tools.ToolChain


class MySFrameBatch(SFrame):

    def configure(self):
        self.xml_doctype = self.xml_doctype + """
<!--
   <ConfigParse NEventsBreak="100000" FileSplit="0" AutoResubmit="0" />
   <ConfigSGE RAM ="2" DISK ="2" Mail="heiner@cern.de" Notification="as" Workdir="workdir"/>
-->
"""
        if os.path.exists(self.cwd + 'workdir'):
            opt = ' -rl --exitOnQuestion'
        else:
            opt = ' -sl --exitOnQuestion'

        self.exe = 'sframe_batch.py' + opt


sframe_tools = ToolChain(
    'SFrameUncerts',
    list(
        SFrame(
            cfg_filename=sys.argv[1],
            xml_tree_callback=set_uncert_func(uncert),
            name='SFrame_' + uncert,
            halt_on_exception=False,
        ) 
        for uncert in sys_uncerts
    )
)


if __name__ == '__main__':
    varial.tools.Runner(sframe_tools)
