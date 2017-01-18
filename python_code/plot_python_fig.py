if __name__ == "__main__":
  import matplotlib.pyplot as plt
  import scipy.io
  import numpy as np

  materials_list = ['espuma_preta_aaa', 'fibra_amarela_aaa', 'la32_aaa', 'la64_aaa', 'melamina_aaa', 'melamina_30amostras2_aaa']

  limit = len(materials_list) - 2
  for material in xrange(0, limit): 
    material_mat = scipy.io.loadmat(materials_list[material] + '.mat')

    ydata = material_mat['ydata']
    ydata = np.ravel(ydata)
    xdata = material_mat['xdata']
    xdata = np.ravel(xdata)

    plt.figure(material)
    for i in xrange(0,12):
      x = xdata[i]
      x = x[0]
      y = ydata[i]
      y = y[0]
      if i <= 6 and i > 2:
        plt.semilogy(x,y, color='red')
      if i >= 7 and i > 2:
        plt.semilogy(x,y, color='black')
      if i == 0:
        plt.semilogy(x,y, 'o', markersize=7, markeredgewidth=1, markeredgecolor='red', markerfacecolor='None')
      if i == 1:
        plt.semilogy(x,y, 'o', markersize=7, markeredgewidth=1, markeredgecolor='black', markerfacecolor='None')

      
    plt.xlabel('frequency')
    plt.ylabel('alpha')
    plt.grid()
    plt.savefig(materials_list[material] + '.png', format='png', dpi=1000)
    #plt.show()

  # fazendo um caso especial para a melamina
  material = 4
  material_mat = scipy.io.loadmat(materials_list[material] + '.mat')

  ydata = material_mat['ydata']
  ydata = np.ravel(ydata)
  xdata = material_mat['xdata']
  xdata = np.ravel(xdata)

  plt.figure(material)
  for i in xrange(0,12):
    x = xdata[i]
    x = x[0]
    y = ydata[i]
    y = y[0]
    if i <= 6 and i > 2:
      plt.semilogy(x,y, color='black')
    if i >= 7 and i > 2:
      plt.semilogy(x,y, color='red')
    if i == 0:
      plt.semilogy(x,y, 'o', markersize=7, markeredgewidth=1, markeredgecolor='red', markerfacecolor='None')
    if i == 1:
      plt.semilogy(x,y, 'o', markersize=7, markeredgewidth=1, markeredgecolor='black', markerfacecolor='None')

    
  plt.xlabel('frequency')
  plt.ylabel('alpha')
  plt.grid()
  plt.savefig(materials_list[material] + '.png', format='png', dpi=1000)
  #plt.show()

  # fazendo um caso especial para a melamina de 30 amostras
  material = 5
  material_mat = scipy.io.loadmat(materials_list[material] + '.mat')

  ydata = material_mat['ydata']
  ydata = np.ravel(ydata)
  xdata = material_mat['xdata']
  xdata = np.ravel(xdata)

  # print ydata[0].shape
  # print ydata[28].shape
  # print ydata[29].shape
  # print ydata[30].shape

  plt.figure(material)
  for i in xrange(0,61):
    x = xdata[i]
    x = x[0]
    y = ydata[i]
    y = y[0]
    if i <= 29 and i > 0:
      plt.semilogy(x,y, 'o', markersize=7, markeredgewidth=1, markeredgecolor='purple', markerfacecolor='None')
    if i >= 30 and i > 0:
      plt.semilogy(x,y, color='black')
    if i == 0:
      plt.semilogy(x,y, 'o', markersize=7, markeredgewidth=1, markeredgecolor='black', markerfacecolor='None')

    
  plt.xlabel('frequency')
  plt.ylabel('alpha')
  plt.grid()
  plt.savefig(materials_list[material] + '.png', format='png', dpi=1000)
  #plt.show()
