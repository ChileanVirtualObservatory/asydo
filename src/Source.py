#SPEED_OF_LIGHT = 299792458.0
#S_FACTOR = 2.3548200450309493820231386529193992754947713787716410 #sqrt(8*ln2) #TODO Preguntar Profe que es esto?


class Source:
    """A source of EM Waves"""

    def __init__(self, log, name, alpha, delta):
        log.write('Source \'' + name + '\' added \n')
        self.log = log
        self.alpha = alpha
        self.delta = delta
        self.name = name
        self.comp = list()

    def addComponent(self, model):
        code = self.name + '-' + str(len(self.comp) + 1)
        self.comp.append(model)
        log.write('Component added to \'' + self.name + '\' (' + code + ')\n')
        log.write('ModelInfo: '+ model.info() + '\n')
# Code passed to Component
#        intens = dict()
#        for mol in mol_list.split(','):
#            abun = random.uniform(conf.base_abun[0], conf.base_abun[1])
#            if mol in ('COv=0', '13COv=0', 'C18O', 'C17O', '13C18O'):
#                abun += conf.base_CO
#            for iso in conf.iso_abun:
#                if iso in mol:
#                    abun *= conf.iso_abun[iso]
#            intens[mol] = abun


# Code moved to Surface
#    def genSurface(self, form, cube):
#        "Create a gaussian surface over a mesh created by x and y axes"
#        sx = form[1]
#        sy = form[2]
#        theta = form[3]
#        r = 3 * math.sqrt(sx ** 2 + sy ** 2)
#        alpha_mesh, delta_mesh = np.meshgrid(cube.alpha_axis, cube.delta_axis, sparse = False, indexing = 'xy')
#        Xc = alpha_mesh.flatten() - self.alpha * np.ones(len(cube.alpha_axis) * len(cube.delta_axis))
#        Yc = delta_mesh.flatten() - self.delta * np.ones(len(cube.alpha_axis) * len(cube.delta_axis))
#        XX = (Xc) * math.cos(theta) - (Yc) * math.sin(theta);
#        YY = (Xc) * math.sin(theta) + (Yc) * math.cos(theta);
#        u = (XX / sx) ** 2 + (YY / sy) ** 2;
#        sol = sx * sy * np.exp(-u / 2) / (2 * math.pi);
#        res = np.transpose(np.reshape(sol, (len(cube.delta_axis), len(cube.delta_axis))))
#        res = res / res.max()
#        return res

    def project(self, cube):
        #self.log.write('Loading visible lines in band ' + cube.band + ' (redshift=' + str(self.rad_vel) + ')\n')
        #db = lite.connect('db/lines.db')
        for component in self.comp:
            log.write(' --> Component name: ' + struct.code + '\n')
            component.project(cube);
            #Moved to Component
            #fwhm = struct.spe_form[1]
            #shape = struct.spe_form[2]
            #freq_init_corr = cube.freq_border[0] * 1000.0 / (1 + self.rad_vel * 1000.0 / SPEED_OF_LIGHT)
            #freq_end_corr = cube.freq_border[1] * 1000.0 / (1 + self.rad_vel * 1000.0 / SPEED_OF_LIGHT)

            #for mol in struct.intens:
            #    log.write("SQL SENTENCE:\n")
            #    select = "SELECT * FROM Lines WHERE SPECIES like '" + mol + "' AND FREQ > " + str(
            #        freq_init_corr) + " AND FREQ < " + str(freq_end_corr)
            #    log.write(select + '\n')
            #    resp = db.execute(select)
            #    linlist = resp.fetchall()
            #    rinte = inten_values[0]
            #    for j in range(len(inten_group)):
            #        if mol in inten_group[j]:
            #            rinte = inten_values[j]
            #    rinte = random.uniform(rinte[0], rinte[1])
            #    for lin in linlist:
            #        trans_temp = lin[5]
            #        temp = np.exp(-abs(trans_temp - self.temp) / self.temp) * rinte
            #        freq = (1 + self.rad_vel * 1000.0 / SPEED_OF_LIGHT) * lin[3]
            #        log.write('E: ' + str(lin[2]) + ' (' + str(lin[1]) + ') at ' + str(freq) + ' Mhz, T = exp(-|' \
            #                  + str(trans_temp) + '-' + str(self.temp) + '|/' + str(self.temp) + ')*' + str(
            #            rinte) + ' = ' + str(temp) + ' K  ')
            #        window = cube.freqWindow(freq / 1000.0, fwhm)
            #        struct.addTransition(mol, str(lin[2]), str(lin[3]), freq, fwhm, temp)
            #        log.write('[W:' + str(window) + ']\n')
            #        sigma = fwhm / S_FACTOR
            #        distro = list()
            #        for idx in range(window[0], window[1]):
            #            distro.append(
            #                np.exp((-0.5 * (cube.freq_axis[idx] - freq / 1000.0) ** 2) / (sigma ** (2 * shape))))
            #        distro = distro / max(distro)
            #        for idx in range(window[0], window[1]):
            #            cube.data[idx] = cube.data[idx] + struct.template * temp * distro[idx - window[0]]

