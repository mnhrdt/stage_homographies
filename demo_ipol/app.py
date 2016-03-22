"""
just another ipol demo
"""

from lib import base_app, build, http, image
from lib.misc import app_expose, ctime
from lib.base_app import init_app
import cherrypy
from cherrypy import TimeoutError
import os.path
import shutil
import time
import math

class app(base_app):
    """demo main class."""

    # IPOL demo system configuration
    title = "Projective Resampling"
    xlink_article = ''

    input_nb = 1          # number of input images
    input_max_pixels = 600*600 # max size (in pixels) of an input image
    input_max_weight = 3 * input_max_pixels  # max size (in bytes)
                                             # of an input file
    input_dtype = '3x8i'  # input image expected data type
    input_ext = '.png'    # input image expected extension (i.e. file format)
    is_test = False       # switch to False for deployment

    def __init__(self):
        """Set up application."""

        # setup the parent class
        base_dir = os.path.dirname(os.path.abspath(__file__))
        base_app.__init__(self, base_dir)

        # select the base_app steps to expose
        # index() is generic
        app_expose(base_app.index)
        app_expose(base_app.input_select)
        app_expose(base_app.input_upload)
        # params() is modified from the template
        app_expose(base_app.params)
        # run() and result() must be defined here

    def build(self):
        return
        #""""Download and compile LSD program."""

        ## store common file path in variables
        #lsd_url = "http://www.ipol.im/pub/art/2012/gjmr-lsd"
        #tgz_name = "lsd_1.6.zip"
        #tgz_file = self.dl_dir + tgz_name
        #build_dir = self.src_dir + "lsd_1.6"
        #bin_file = build_dir + "/lsd"
        #prog_file = self.bin_dir + "lsd"
        #log_file = self.base_dir + "build.log"

        ## get the latest source archive
        #build.download(lsd_url+"/"+tgz_name, tgz_file)

        ## test if the dest file is missing, or too old
        #if (os.path.isfile(prog_file)
        #    and ctime(tgz_file) < ctime(prog_file)):
        #    cherrypy.log("not rebuild needed",
        #                 context='BUILD', traceback=False)
        #else:
        #    # extract the archive
        #    build.extract(tgz_file, self.src_dir)

        #    # build the program
        #    build.run("make -C %s" % build_dir, stdout=log_file)

        #    # save into bin dir
        #    if os.path.isdir(self.bin_dir):
        #        shutil.rmtree(self.bin_dir)
        #    os.mkdir(self.bin_dir)
        #    shutil.copy(bin_file, prog_file)

        #    # cleanup the source dir
        #    shutil.rmtree(self.src_dir)
        #return

    def set_default_gui_state(self):
        img = image(self.work_dir + 'input_0.png')
        w = img.size[0]
        h = img.size[1]
        m = 50
        ow = w + 2 * m
        oh = h + 2 * m

        self.cfg['param']['m'] = m
        self.cfg['param']['w'] = w
        self.cfg['param']['h'] = h
        self.cfg['param']['ow'] = ow
        self.cfg['param']['oh'] = oh
        self.cfg['param']['c'] = str([[(0+m+int(w*0.4)),0+m+int(h*0.4)],[w-1+m,0+m],[w-1+m,h-1+m],[0+m,h-1+m]])
        self.cfg['param']['p'] = str([[0,0],[w-1,0],[w-1,h-1],[0,h-1]])
        self.cfg['param']['hit'] = -1
        self.cfg['param']['state'] = "plain"
        return

    @cherrypy.expose
    @init_app
    def params(self, newrun=False, msg=None):
        """Parameter handling"""

        # if a new experiment on the same image, clone data
        if newrun:
            self.clone_input()

        self.set_default_gui_state()
        self.cfg['param']['has_gallery'] = "gallery_old"
        return self.tmpl_out('params.html')

    @cherrypy.expose
    @init_app
    def wait(self, **kwargs):
        """Input handling and run redirection."""

        # otherwise, run the algorithm
        http.refresh(self.base_url + 'run?key=%s' % self.key)
        return self.tmpl_out("wait.html")

    @cherrypy.expose
    @init_app
    def run(self):
        """Run the algorithm."""

        try:
            run_time = time.time()
            self.run_algo()
            self.cfg['info']['run_time'] = time.time() - run_time
        except TimeoutError:
            return self.error(errcode='timeout',
                              errmsg="Try again with simpler images.")
        except RuntimeError:
            return self.error(errcode='runtime',
                              errmsg="Something went wrong with the program.")
        http.redir_303(self.base_url + 'result?key=%s' % self.key)

        # archive
        if self.cfg['meta']['original']:
            ar = self.make_archive()
            ar.add_file("input_0.orig.png","uploaded_image.png")
            ar.add_info({"run time (s)" : self.cfg['info']['run_time']})
            ar.save()

        return self.tmpl_out("run.html")

    def update_alphadots(self):
        ow = self.cfg['param']['ow']
        oh = self.cfg['param']['oh']
        c = eval(self.cfg['param']['c'])
        hit = int(self.cfg['param']['hit'])
        colors = ["gray", "gray", "gray", "gray"]
        if (hit >= 0 and hit < 4):
            colors[hit] = "red"
        cmdline = ["alphadots", "-p", str(ow), str(oh), "dots.png",
                           str(c[0][0]), str(c[0][1]), colors[0],
                           str(c[1][0]), str(c[1][1]), colors[1],
                           str(c[2][0]), str(c[2][1]), colors[2],
                           str(c[3][0]), str(c[3][1]), colors[3]
                  ]
        #self.wait_proc(self.run_proc(["/bin/echo", str(cmdline)]))
        self.wait_proc(self.run_proc(cmdline))
        return

    def update_one_view(self,order,outname):
        ow = self.cfg['param']['ow']
        oh = self.cfg['param']['oh']
        c = eval(self.cfg['param']['c'])
        p = eval(self.cfg['param']['p'])

        if (order == 1):
            cmd = "warp_decomp_ripmap"
            cmdline = [cmd, "input_0.png",
                       str(c[0][0]), str(c[0][1]),
                       str(c[1][0]), str(c[1][1]),
                       str(c[2][0]), str(c[2][1]),
                       str(c[3][0]), str(c[3][1]),
                       str(p[0][0]), str(p[0][1]),
                       str(p[1][0]), str(p[1][1]),
                       str(p[2][0]), str(p[2][1]),
                       str(p[3][0]), str(p[3][1]),
                       str(ow), str(oh)
                      ]
        else:
            cmd = "viho_demo"
            cmdline = [cmd, "-m", str(order),
                       "input_0.png", str(outname), str(ow), str(oh),
                       str(c[0][0]), str(c[0][1]),
                       str(c[1][0]), str(c[1][1]),
                       str(c[2][0]), str(c[2][1]),
                       str(c[3][0]), str(c[3][1]),
                       str(p[0][0]), str(p[0][1]),
                       str(p[1][0]), str(p[1][1]),
                       str(p[2][0]), str(p[2][1]),
                       str(p[3][0]), str(p[3][1])
                      ]
        self.wait_proc(self.run_proc(["/bin/echo", str(cmdline)]))
        self.wait_proc(self.run_proc(cmdline))

    def update_canvas(self):
        self.update_one_view(0, "o_nn.png")
        return

    def run_algo(self):
        """Core algorithm runner, it could also be called by a batch processor,
           this one needs no parameter.
        """

        # dots
        self.update_alphadots()

        # canvas
        self.update_canvas()

        # run the program
        #self.wait_proc(self.run_proc(['do_everything.sh']))
        self.run_comparison()
        return

    def did_hit(self, x, y, i):
        c = eval(self.cfg['param']['c'])
        ix = int(c[i][0])
        iy = int(c[i][1])
        return math.hypot(x - ix, y - iy) < 20

    def compute_hit(self, x, y):
        if self.did_hit(x, y, 0): return 0
        if self.did_hit(x, y, 1): return 1
        if self.did_hit(x, y, 2): return 2
        if self.did_hit(x, y, 3): return 3
        return -1

    @cherrypy.expose
    @init_app
    def reset(self, **kwargs):
        self.set_default_gui_state()
        self.update_alphadots()
        self.update_canvas()
        self.cfg['param']['has_gallery'] = "gallery_old"
        return self.tmpl_out("result.html", height=1000)

    @cherrypy.expose
    @init_app
    def click(self, **kwargs):
        print("\nIN CLICK\n\n")
        print(kwargs)
        px = int(kwargs['canvas.x'])
        py = int(kwargs['canvas.y'])
        state = self.cfg['param']['state']
        if state == "plain":
            hit = self.compute_hit(px,py)
            self.cfg['param']['px'] = px
            self.cfg['param']['py'] = px
            self.cfg['param']['hit'] = hit
            print("hit = %d\n" % hit)
            if hit >= 0:
                self.update_alphadots()
                self.cfg['param']['state'] = "wait_other"
        elif state == "wait_other":
            print("new pos = %d %d\n" % (px, py))
            hit = int(self.cfg['param']['hit'])
            if hit >= 0:
                c = eval(self.cfg['param']['c'])
                c[hit] = [px, py]
                self.cfg['param']['c'] = str(c)
                self.cfg['param']['state'] = "plain"
                self.cfg['param']['hit'] = "-1"
                self.update_alphadots()
                self.update_canvas()
                self.cfg['param']['has_gallery'] = "gallery_old"
        return self.tmpl_out("result.html", height=1000)

    def run_comparison(self):
        self.update_alphadots()
        self.update_one_view(0, "out_nn.png")
        self.update_one_view(-1, "out_ndots.png")
        self.update_one_view(2, "out_bil.png")
        self.update_one_view(3, "out_bic.png")

        run_time = time.time()
        self.update_one_view(1, "out_final.png")
        self.cfg['info']['run_time'] = time.time() - run_time

        self.cfg['param']['has_gallery'] = "gallery_new"
        return

    @cherrypy.expose
    @init_app
    def rerun(self, **kwargs):
        self.run_comparison()
        return self.tmpl_out("result.html", height=1000)

    @cherrypy.expose
    @init_app
    def result(self):
        """Display the algorithm result."""

	h = self.cfg['param']['h']

        return self.tmpl_out("result.html", height=1000)
