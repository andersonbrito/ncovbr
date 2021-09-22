# -*- coding: utf-8 -*-

import pycountry_convert as pyCountry
import pycountry
import pandas as pd
import argparse
from uszipcode import SearchEngine


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Reformat metadata file by adding column with subcontinental regions based on the UN geo-scheme",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", required=True, help="Nextstrain metadata file")
    parser.add_argument("--geoscheme", required=True, help="XML file with geographic classifications")
    parser.add_argument("--output", required=True, help="Updated metadata file")
    args = parser.parse_args()

    metadata = args.metadata
    geoscheme = args.geoscheme
    output = args.output

    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/ITpS/projetos_itps/dashboard/nextstrain/run5_20210920_template/'
    # metadata = path + 'pre-analyses/metadata_filtered.tsv'
    # geoscheme = path + 'config/geoscheme.tsv'
    # output = path + 'pre-analyses/metadata_geo.tsv'


    # get ISO alpha3 country codes
    isos = {}
    def get_iso(country):
        global isos
        if country not in isos.keys():
            try:
                isoCode = pyCountry.country_name_to_country_alpha3(country, cn_name_format="default")
                isos[country] = isoCode
            except:
                try:
                    isoCode = pycountry.countries.search_fuzzy(country)[0].alpha_3
                    isos[country] = isoCode
                except:
                    isos[country] = ''
        return isos[country]

    # parse subcontinental regions in geoscheme
    scheme_list = open(geoscheme, "r").readlines()[1:]
    geoLevels = {}
    c = 0
    for line in scheme_list:
        if not line.startswith('\n'):
            id = line.split('\t')[2]
            type = line.split('\t')[0]
            members = line.split('\t')[5].split(',')  # elements inside the subarea

            if type == 'region':
                for country in members:
                    iso = get_iso(country.strip())
                    geoLevels[iso] = id

            # parse subnational regions for countries in geoscheme
            if type == 'country':
                for state in members:
                    if state.strip() not in geoLevels.keys():
                        geoLevels[state.strip()] = id

            for elem in members:
                if elem.strip() not in geoLevels.keys():
                    geoLevels[elem.strip()] = id

    # open metadata file as dataframe
    dfN = pd.read_csv(metadata, encoding='utf-8', sep='\t')
    try:
        dfN.insert(4, 'region', '')
    except:
        pass
    dfN['region'] = dfN['code'].map(geoLevels) # add 'column' region in metadata
    dfN['br_region'] = ''

    # convert sets of locations into sub-locations
    print('\nApplying geo-schemes...')
    dfN.fillna('', inplace=True)

    for idx, row in dfN.iterrows():
        country = dfN.loc[idx, 'country']

        # assign BR region
        if country not in ['Brazil']:
            dfN.loc[idx, 'br_region'] = 'Other regions'

        if country == 'Brazil' and dfN.loc[idx, 'br_region'] == '':
            dfN.loc[idx, 'br_region'] = dfN.loc[idx, 'division']

        # divide country into subnational regions
        division = dfN.loc[idx, 'division']
        if division not in ['', 'unknown']:
            if division in geoLevels.keys():
                dfN.loc[idx, 'country'] = geoLevels[dfN.loc[idx, 'division']]

        print('Processing metadata for... ' + row['strain'])

        # country = dfN.loc[idx, 'country']
        # division = dfN.loc[idx, 'division']
        #
        # if country == 'Brazil':
        #     if division not in ['', None]:
        #         dfN.loc[idx, 'country'] = custom_geolevels[division]
        #         dfN.loc[idx, 'br_region'] = division
        #     else:
        #         division = 'Unknown'
        #         dfN.loc[idx, 'country'] =custom_geolevels[division]
        # else:
        #     dfN.loc[idx, 'br_region'] = 'Other regions'


    dfN = dfN.drop_duplicates(subset=['strain'])
    dfN.to_csv(output, sep='\t', index=False)

print('\nMetadata file successfully reformatted applying geo-scheme!\n')
