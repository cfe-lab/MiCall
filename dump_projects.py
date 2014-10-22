import json
import requests

import settings

def main():
    dump = {}
    with requests.Session() as session:
        response = session.post(settings.qai_path + "/account/login",
                                data={'user_login': settings.qai_user,
                                      'user_password': settings.qai_password})
        # Would be nice to check response code, but it's always 200
        if 'loginform' in response.text:
            exit('Login failed, check qai_user in settings.py')
        
        regions = session.get(settings.qai_path + "/lab_miseq_regions.json?mode=dump",
                              auth=(settings.qai_user,
                                    settings.qai_password))
        projects = session.get(settings.qai_path + "/lab_miseq_projects.json?mode=dump",
                              auth=(settings.qai_user,
                                    settings.qai_password))
        dump['regions'] = regions.json()
        dump['projects'] = projects.json()
        
        
    with open("projects.json", "w") as f:
        json.dump(dump, f, sort_keys=True, indent=2, separators=(',', ': '))
    
    print "Done."

if __name__ == "__main__":
    main()
