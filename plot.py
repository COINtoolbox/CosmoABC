import pylab as plt
import kdestats as kde

########################################
#root = 'SMC_ABC_SZ_'
root = "SimpleTest_"

########################################


for i in range( 18, 19 ): 

    path = root + str( i ) + '.dat'

    op1 = open( path, 'r')
    lin1 = op1.readlines()
    op1.close()

    data1 = [ elem.split() for elem in lin1 ]

    om = [ float( item[0] ) for item in data1[1:] ]
   # w = [ float( item[1] ) for item in data1[1:] ]
    sig8 = [ float( item[1] ) for item in data1[1:] ]

    plt.figure(figsize=(10,10))
  #  plt.subplot(3,2,1)
  #  kdehistA = kde.kdehist2(om, w, [30, 30])
  #  clevelsA = kde.confmap(kdehistA[0], [.6827,.9545,.9973])
  #  plt.contour(kdehistA[1], kdehistA[2], kdehistA[0], clevelsA)
  #  plt.xlabel( 'om' )
  #  plt.ylabel( 'w' )
  #  plt.xlim(0,1)
  #  plt.ylim( -3.0,1.0)
  #  plt.title('Rejection loop : ' + str( i ) )

    plt.subplot( 3,1,1)
    plt.hist( om )
    plt.xlabel( 'om' )
#    plt.xlim( 0, 1)
    plt.ylabel( 'number count' )

    plt.subplot( 3,1,3)
    kdehistB = kde.kdehist2(om, sig8, [30, 30])
    clevelsB = kde.confmap(kdehistB[0], [.6827,.9545,.9973])
    plt.contour(kdehistB[1], kdehistB[2], kdehistB[0], clevelsB)
    plt.xlabel( 'om' )
    plt.ylabel( 'sig8' )
   # plt.xlim(0,1)
   # plt.ylim(0.0,1)

  #  plt.subplot(3,2,4)
  #  plt.hist( w )
  #  plt.xlabel( 'w' )
  #  plt.xlim( -3.0,1.0)
  #  plt.ylabel( 'number count' )

  #  plt.subplot( 3,2,5)
  #  kdehistC = kde.kdehist2( w, sig8, [30, 30])
  #  clevelsC = kde.confmap(kdehistC[0], [.6827,.9545,.9973])
  #  plt.contour(kdehistC[1], kdehistC[2], kdehistC[0], clevelsC)
  #  plt.xlabel( 'w' )
  #  plt.ylabel( 'sig8' )
  #  plt.xlim( -2.0,1.0)
  #  plt.ylim(0.6,1.0)

    plt.subplot(3,1,2)
    plt.hist( sig8 )
    plt.xlabel( 'sig8' )
   # plt.xlim( 0.0,1.0)
    plt.ylabel( 'number count' )
    
    plt.savefig(root + str( i )  )


#plt.show()
plt.close()
