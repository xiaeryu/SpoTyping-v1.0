import urllib
import urllib2

class Crawler(object):
    DEBUG = False
    TIMEOUT = 500

    def __init__(self, user=None):
        #enable cookie
        self.opener = urllib2.build_opener(urllib2.HTTPCookieProcessor())
        pass

    def get(self, url):
        response = self.opener.open(url, timeout=self.TIMEOUT)
        return response.read()

    def post(self, url, data):
        req = urllib2.Request(url)
        data = urllib.urlencode(data)
        response = self.opener.open(req, data, timeout=self.TIMEOUT)
        return response.read()

    @staticmethod
    def join_str(iterable=[], seperator="," ):
        new_items = []
        for item in iterable:
            if isinstance(item, basestring):
                new_items.append(item)
        return seperator.join(new_items)

    # def load(self, url, filename):
    #     if not self.DEBUG:
    #         response = self.get(url)
    #         with open(filename, "w") as file:
    #             file.write(response)
    #     else:
    #         with open(filename) as file:
    #             response = file.readlines()
    #     return response

    def run(self):
        pass
