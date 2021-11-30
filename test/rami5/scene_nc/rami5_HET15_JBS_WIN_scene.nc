CDF       
      column        layer      >   layer_interface    ?   
wavelength           	   title         EForest description for the RAMI5 J�rvselja Birch Stand (Winter) scene      scene_id      HET15_JBS_WIN      scene         J�rvselja Birch Stand (Winter)     summary      VThis file contains a description of a RAMI5 scene ready to be read in by the SPARTACUS-Surface
radiative transfer model. In order that this model can treat tree trunks, we treat these forests as
vegetated urban areas, where the trunks are treated as buildings. The radiative effects of other
woody material (branches) has been added to the vegetation properties. Changes in version 0.8:
crown radius is 0.9x max leaf distance from tree axis; crowns < 2m from the ground are projected to
the ground to prevent erroneous horizontal transport on reflection; veg_fsd is provided as an
explicit profile.     institution       2European Centre for Medium-Range Weather Forecasts     creator_name      Robin Hogan    creator_email         r.j.hogan@ecmwf.int    creation_date         30-Nov-2021    software_version      0.8          
wavelength                 	long_name         
Wavelength     units         nm     comment       UA dummy wavelength entry of zero is added at the end, for which all albedos are zero.         8  �   surface_type                	long_name         Surface type   
definition        ?0: Flat
1: Forest
2: Urban
3: Woody forest (or vegetated urban)         0   cos_solar_zenith_angle                  	long_name          Cosine of the solar zenith angle   units         1           4   top_flux_dn_sw                     	long_name         Top-of-canopy downwelling flux     units         W m-2         8  8   top_flux_dn_direct_sw                      	long_name         %Top-of-canopy downwelling direct flux      units         W m-2         8  p   nlayer                  	long_name         Number of active layers         �   height                     	long_name         Height of layer interfaces     units         m         �  �   building_fraction                      	long_name         Trunk cover fraction   units         1      comment       dThis variable is computed analytically from the number of trees falling withing the analysis domain.      �  �   veg_fraction                   	long_name         Vegetation fraction    units         1      comment       �This variable combines tree crowns and bare branches, and is computed from a rasterization of the
analysis domain in which circular trees are placed.         �  �   building_scale                     	long_name         Trunk effective diameter   units         m      comment       �This is the diameter of tree trunks in an equivalent forest in which all tree trunks have the same
diameter, and have the same areal fraction and surface area as the trees in the real forest.       �  �   	veg_scale                      	long_name         Vegetation horizontal scale    units         m      comment       �This is the variable S used in Eq. 20 of Hogan et al. (GMD 2018) relating the perimeter length of
the vegetation and the areal fraction.      �  �   veg_extinction                     	long_name         !Vegetation extinction coefficient      units         m-1    comment       �This variable is the wavelength-independent extinction coefficient of the vegetated fraction of
each layer, combining the contribution from foliage and tree branches (but not trunks) assuming the
surfaces of both to be randomly oriented.         �  �   veg_contact_fraction                   	long_name         Vegetation contact fraction    units         1      comment       |This is the fraction of the perimeter of the trunks in each layer that is in contact with
vegetation, rather than clear air.      �  �   ground_sw_albedo                   	long_name         Ground albedo      units         1         8   x   wall_sw_albedo                        	long_name         Trunk albedo   units         1        �   �   roof_sw_albedo                        	long_name         Trunk narrowing albedo     units         1      comment      BThe area coverage of tree trunks reduces with height as the trunks get narrower. This reduction is
represented in a stepwise fashion, with a constant trunk diameter within individual layers. This
variable represents the albedo of the small horizontal surface exposed between two layers with
different trunk cover fraction.       �  .@   
veg_sw_ssa                        	long_name         #Vegetation single scattering albedo    units         1      comment       4Scattering by vegetation is assumed to be isotropic.     �  ;�   foliage_extinction                     	long_name         Foliage extinction coefficient     units         m-1    comment      wThis variable is the wavelength-independent extinction coefficient of the vegetated fraction of
each layer, but containing only the contribution from foliage and not branches or trunks. This is
not used by SPARTACUS-Surface, but can be used (with foliage_sw_ssa) after running
SPARTACUS-Surface to work out the amount of radiation absorbed by foliage and branches separately.       �  I`   foliage_sw_ssa                        	long_name          Foliage single scattering albedo   units         1      comment      This variable is like veg_sw_ssa but for the foliage only (no branches). It is not used by
SPARTACUS-Surface, but can be used (with foliage_extinction) after running SPARTACUS-Surface to
work out the amount of radiation absorbed by foliage and branches separately.     �  JX   veg_fsd                    	long_name         3Vegetation extinction fractional standard deviation    units         1      comment      �This variable is the fractional standard deviation (FSD) of vegetation extinction within the crowns
of a tree canopy, representing the clumpiness of leaf and branch distributions. It has been
calculated here as an average of the assumed foliage FSD of 0.7 and branch FSD of 1.4 weighted
according to the contribution of the foliage and branch extinction to the total vegetation
extinction.         �  W�C�y�C�9�D�D&P D*c3D1FfD<��DP` DXffD�P D�P E� E	p      �?   ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�   >�    ?   ?�  ?�  @   @   @@  @`  @�  @�  @�  @�  @�  @�  @�  @�  A   A  A  A  A   A(  A0  A8  A@  AH  AP  AX  A`  Ah  Ap  Ax  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  ;[�M;ND�;M�;LT~;Kþ;J�;I��;I;H;G�;F;E��;D�;D	Q;C8X;Bj;>��;;0*;6U;/�};+F';'L;#a];X;S; �;��;	ZB;6�:���:�]:���:��:��0:Ē�:���:��@:�-�:Vw�:+��:��9�Y�9�89i�t9)&8��8��(8_Iv8!��7��17�L�7��'71��6�_56��6��5`��4n��2�{�0 �@$��d    <�֑<�֑=#=#=#=#=#=#='k=��=��=O�v=;X�=��b=�ki>\�>b�>_"
>X�>]�T>s#�>jU4>^��>UE�>kx�>��<>��.>�%\>�7�>�B�>�^>���>ќ�>�	�?�e?0�?�?��?
�t?�s>��>�G�>�8�>�.	>�X@>�.>��>�X.>��>��">�.2>�|>ظ>��=��=ڸv=�
<��<���<���<<��    >a�M>Z�u>Y�E>Y��>Y5�>X�m>X{�>XU�>X?�>X!�>W�>W�!>W`�>WL3>WE�>W4�>U�^>T�R>SX4>Q�>O�M>O[c>Nh�>M�>M>M��>K,>I�>F�y>E�>D�A>E.C>F�,>F�>D6><��>2�6>&�2>��>
0�=���=��
=��[=���=��8=�|9=��*=|7H=�@=���=�}0=�Y-=t�=`Њ=:�7=��<�֢<,pJ;s�a:o�/�W\    @",@",'@X̷@X̽@X��@X��@X��@X��@U,@E=�@)�@7}@5Ȳ@F"�@9�5@h�n@mAD@k��@a��@\n+@]y�@_
�@\�N@R�z@:R�@F��@G@V@C�<@@�Z@I�@<�Y@C^�@+�@.o|@4U�@9�8@<��@8��@6�@52@9��@9�V@7 �@7!�@78�@D $@Cv�@B��@>�e@;�@4l+@WW�@c.@`c@W�%@Vue@x�S@g�@T��@?�@�G    4X�4[��9��`<A��<��S<�%=��<�@j=�`.=��>z�=$!�=N(=���=�8�<��<�œ=e�=��=~Ke=Mc�=I=N%�=v��=�*�=���=��=�n=�h�=�i�=�D�=�C�=���=�|K=���=��=��'=��=�)9=�b=�=�E|=۟�=�o�=ͻg=��=��-=��[=o�=<x<= �=��=��%=Y�=9��=m2=�Q�=��=S�C<ꍁ;n�    <��><�?<��)<�75<�"<�&�<˥�<�l�<��E=<?=2}U=/�=^E�=�c>h
>c��>Z�,>O��>nt>`��>R�)>B��>�ߍ>�i�>ğ�>�΃>�A^>��/>�1�>焖?��>�'?tw9?x�t?x\?x}S?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�      ?a�a?e١?h1j?i�?i�<?g'$?f�?[MI?Y�>��=9�h<���=Z(5    =Ǻ
=�>R�>��>��>tR�>��y>Û�>��y?o.>��>ϵ9>�C    =�=��>�>1v>>sѡ>�_�>�?6>�{6?i\>�:�>�=�>��p    =��|=�ט>�	>�>�U>s�<>�m:>�P(>Ό~?h�>�&�>�1V>��H    =��1=���>�>�+>��>s��>�s�>�Xu>Δ�?h�>�>�,�>���    =�æ=Գ�>��>ܸ>�e>s�*>�y�>�_�>Μ)?ig>�>�*�>���    =š$=ԉF>��>��>��>s�>���>�n�>Ϋ�?i]>�	>�!>Ͻ�    =�zC=�XA>��>�>c�>s� >���>ĂN>οi?i�>���>��>ϴ�    =�L�=��>|�>N�>+�>s��>���>Ė�>��?i�>��x>��>ϩ�    =�4<=���>_�>�>��>s�p>���>į<>���?l�>��D>�6>ϥ�    =�H=���>?O>��>�}>s��>��[>��>��?oq>��>�^>Ϡn    =��=Ӹe>&�>��>�>szM>��>���>��?r�>��B>��>Ϟ�    =� =Ӭ>K>��>{�>sv�>���>��>�)�?t�>��8>�u>Ϟ7    =��=ӗ4>
�>�o>^�>sqd>��L>��5>�7m?v>���>� .>ϛ�    =��=ӆE>�0>j<>7�>slf>�3>�>�MV?y�>���>�.>Ϝ    =��V=�_>أ>8,>S>sb>� �>�'S>�e�?|g>��;>��u>ϗ    =ĕ=��>��>�c>�a>sO0>�<>�H�>χ�?~>�w�>��>ϊ�    =��=�F)>$\>+>�>s�>�t+>ŏ�>���?}5>�*�>�­>�a�    =�S=с�>�>`}>�>r�>���>���>�%2?��>��O>Ν�>�>�    =¾f=н>�>�->:�>r��>��M>�2�>�u ?>��U>�p�>��    =�\=�ף> �>��>HY>rgX>�7`>��>��N?u�>��>�0>�֍    =�2�=λ�=���>�H>0�>r�>��b>��h>�?nz>���>��>Γi    =���=�@<=��f>Fs>��>q�>��B>�W/>ѝ	?j4>��>͐0>�C�    =���=��=�
<>p>{e>qZ�>�:�>���>��?d�>�>�<+>���    =��G=�nR=�G!>��>$�>p��>���>�:�>҄�?_>���>��>ͨl    =�r�=��[=�j
>
b>�d>p��>��>ȵ5>�F?Yt>�jQ>̊>�U�    =��p=�$�=�9�>ɰ>�>p%e>�h�>�FT>ӕ?S�>�>� z>��$    =�K=Ĳ=�0s>��>
��>om�>��%>���>�Fh?9h>�>�p�>�T*    =���=��=��>��>	�>nϩ>�C�>�^@>ԯ?>ￆ>���>��%    =���=��%=�I�>�y>��>nG>�m�>ʕ]>���?��>���>�E�>�:    =��=�L/=��}>��>�>m��>��3>��O>�N?��>�*h>ɹ�>ʵT    =��H=��]=�ߛ>@*>7�>m-^>��>�>�S�?��>�P>�"+>�%�    =���=�C�=�כ=���>�~>l�\>��>�e�>ղs?n">�z�>ȓo>ɠ    =��=�E�=�XN=��>��>k��>�c�>��.>�@??�>�j�>�۳>��Z    =��s=�q{=��=��g> >k$�>���>�d>�e;? ��>�Gz>��>�2    =��@=�!�=�R�=�W=�|�>j^>���>�>>�X�? ��>� 0>�7 >�`5    =��=��p=��=��=�J2>i�C>�!�>ˊ~>��? f>��S>�9`>�ef    =�D=�\�=��=�8�=��>h��>�t�>ʶ<>��R>�ٯ>�w>�>�>G    =�v�=�"H=ᣵ=��=�>>g�>� F>�'
>�H>���>�r�>�=�>�l    =��R=��=�W�=�5�=���>g!�>�o�>�t�>ӊ#>�ɬ>�S^>�P�>À    =�L=���=��=���=�F�>fG�>���>ȟV>Ҧ�>�� >��>�EC>�t�    =��2=��9=��=�H�=��>e:�>�� >Ǫ�>ѡ�>���>ၸ>��u>�-�    =��)=��=�a=�c�=��\>c��>�Jy>��R>й
>�FW>߆>�j>��    =��=�/r=�E�=��_=��>b��>���>� e>���>���>ݲ�>���>�-�    =��X=�JV=��`=��=��>aq>��>�Ho>��>�#l>��<>���>�Î    =�9�=�.r=�qp=�=�Ԕ>`d�>��n>���>Θ0>��>�S�>�Q�>���    =��G=�=�/�=�S�=�w<>_��>�$�>�4�>��>��J>�G�>�t�>��z    =�^K=�Y$=�$�=�Ж=��@>^�t>��>���>�m>���>�#>�y@>���    =�(=�sr=�	2=��=��S>]�l>�Ӷ>3>�*I>�)�>�^Q>��>�^�    =���=�U�=ڐX=�{�=� >[3>�ZP>���>��>��>ҡ >��o>�*X    =�0~=�#v=�f�=���=��>V�D>���>��>�&>>�,A>�߻>�M$>���    =���=�U�=�=�=�Of=�r�>Q�>�f�>��_>��>�)>�a�>�8D>���    =��(=�e�=�Z�=��=�E>L��>��7>�30>���>���>�e�>�g�>���    =�5�=��=��=��=�G�>Go�>�j�>��%>��>�<�>�VX>���>�چ    =��=��&=��=�;=�.�>A1Y>�a�>�uy>�!�>��>�_>��>��@    =��Q=�t=�2�=�t�=�1>=�J>���>��(>�A>��>��]>�">�Q�    =��=���=�B�=�-=�K�>=eh>�x&>��h>��>���>�xq>���>��    =��=��=��u=�	�=���>>]9>�x�>�ܼ>�R�>�J�>��Q>�"�>�T�    =Ơa=��=�#=���=�u(>@�z>�Ĕ>��>�S�>�x>�/>��>�+e    =�܏=���=Ւ�=��=��>E�>�Z�>�^>>�Tg>�~�>���>���>�׏    =���=�%�=�:�=�E=紽>H��>�L>�=�>�}T>�� >�_$>�j >���    =�I�=��j=�J=ˑ�=���>pF�>�.s>��>���?	��>��G>ϐ�>�Ib                                                            =Ǻ
=�>R�>��>��>tR�>��y>Û�>��y?o.>��>ϵ9>�C    =�=��>�>1v>>sѡ>�_�>�?6>�{6?i\>�:�>�=�>��p    =��|=�ט>�	>�>�U>s�<>�m:>�P(>Ό~?h�>�&�>�1V>��H    =��1=���>�>�+>��>s��>�s�>�Xu>Δ�?h�>�>�,�>���    =�æ=Գ�>��>ܸ>�e>s�*>�y�>�_�>Μ)?ig>�>�*�>���    =š$=ԉF>��>��>��>s�>���>�n�>Ϋ�?i]>�	>�!>Ͻ�    =�zC=�XA>��>�>c�>s� >���>ĂN>οi?i�>���>��>ϴ�    =�L�=��>|�>N�>+�>s��>���>Ė�>��?i�>��x>��>ϩ�    =�4<=���>_�>�>��>s�p>���>į<>���?l�>��D>�6>ϥ�    =�H=���>?O>��>�}>s��>��[>��>��?oq>��>�^>Ϡn    =��=Ӹe>&�>��>�>szM>��>���>��?r�>��B>��>Ϟ�    =� =Ӭ>K>��>{�>sv�>���>��>�)�?t�>��8>�u>Ϟ7    =��=ӗ4>
�>�o>^�>sqd>��L>��5>�7m?v>���>� .>ϛ�    =��=ӆE>�0>j<>7�>slf>�3>�>�MV?y�>���>�.>Ϝ    =��V=�_>أ>8,>S>sb>� �>�'S>�e�?|g>��;>��u>ϗ    =ĕ=��>��>�c>�a>sO0>�<>�H�>χ�?~>�w�>��>ϊ�    =��=�F)>$\>+>�>s�>�t+>ŏ�>���?}5>�*�>�­>�a�    =�S=с�>�>`}>�>r�>���>���>�%2?��>��O>Ν�>�>�    =¾f=н>�>�->:�>r��>��M>�2�>�u ?>��U>�p�>��    =�\=�ף> �>��>HY>rgX>�7`>��>��N?u�>��>�0>�֍    =�2�=λ�=���>�H>0�>r�>��b>��h>�?nz>���>��>Γi    =���=�@<=��f>Fs>��>q�>��B>�W/>ѝ	?j4>��>͐0>�C�    =���=��=�
<>p>{e>qZ�>�:�>���>��?d�>�>�<+>���    =��G=�nR=�G!>��>$�>p��>���>�:�>҄�?_>���>��>ͨl    =�r�=��[=�j
>
b>�d>p��>��>ȵ5>�F?Yt>�jQ>̊>�U�    =��p=�$�=�9�>ɰ>�>p%e>�h�>�FT>ӕ?S�>�>� z>��$    =�K=Ĳ=�0s>��>
��>om�>��%>���>�Fh?9h>�>�p�>�T*    =���=��=��>��>	�>nϩ>�C�>�^@>ԯ?>ￆ>���>��%    =���=��%=�I�>�y>��>nG>�m�>ʕ]>���?��>���>�E�>�:    =��=�L/=��}>��>�>m��>��3>��O>�N?��>�*h>ɹ�>ʵT    =��H=��]=�ߛ>@*>7�>m-^>��>�>�S�?��>�P>�"+>�%�    =���=�C�=�כ=���>�~>l�\>��>�e�>ղs?n">�z�>ȓo>ɠ    =��=�E�=�XN=��>��>k��>�c�>��.>�@??�>�j�>�۳>��Z    =��s=�q{=��=��g> >k$�>���>�d>�e;? ��>�Gz>��>�2    =��@=�!�=�R�=�W=�|�>j^>���>�>>�X�? ��>� 0>�7 >�`5    =��=��p=��=��=�J2>i�C>�!�>ˊ~>��? f>��S>�9`>�ef    =�D=�\�=��=�8�=��>h��>�t�>ʶ<>��R>�ٯ>�w>�>�>G    =�v�=�"H=ᣵ=��=�>>g�>� F>�'
>�H>���>�r�>�=�>�l    =��R=��=�W�=�5�=���>g!�>�o�>�t�>ӊ#>�ɬ>�S^>�P�>À    =�L=���=��=���=�F�>fG�>���>ȟV>Ҧ�>�� >��>�EC>�t�    =��2=��9=��=�H�=��>e:�>�� >Ǫ�>ѡ�>���>ၸ>��u>�-�    =��)=��=�a=�c�=��\>c��>�Jy>��R>й
>�FW>߆>�j>��    =��=�/r=�E�=��_=��>b��>���>� e>���>���>ݲ�>���>�-�    =��X=�JV=��`=��=��>aq>��>�Ho>��>�#l>��<>���>�Î    =�9�=�.r=�qp=�=�Ԕ>`d�>��n>���>Θ0>��>�S�>�Q�>���    =��G=�=�/�=�S�=�w<>_��>�$�>�4�>��>��J>�G�>�t�>��z    =�^K=�Y$=�$�=�Ж=��@>^�t>��>���>�m>���>�#>�y@>���    =�(=�sr=�	2=��=��S>]�l>�Ӷ>3>�*I>�)�>�^Q>��>�^�    =���=�U�=ڐX=�{�=� >[3>�ZP>���>��>��>ҡ >��o>�*X    =�0~=�#v=�f�=���=��>V�D>���>��>�&>>�,A>�߻>�M$>���    =���=�U�=�=�=�Of=�r�>Q�>�f�>��_>��>�)>�a�>�8D>���    =��(=�e�=�Z�=��=�E>L��>��7>�30>���>���>�e�>�g�>���    =�5�=��=��=��=�G�>Go�>�j�>��%>��>�<�>�VX>���>�چ    =��=��&=��=�;=�.�>A1Y>�a�>�uy>�!�>��>�_>��>��@    =��Q=�t=�2�=�t�=�1>=�J>���>��(>�A>��>��]>�">�Q�    =��=���=�B�=�-=�K�>=eh>�x&>��h>��>���>�xq>���>��    =��=��=��u=�	�=���>>]9>�x�>�ܼ>�R�>�J�>��Q>�"�>�T�    =Ơa=��=�#=���=�u(>@�z>�Ĕ>��>�S�>�x>�/>��>�+e    =�܏=���=Ւ�=��=��>E�>�Z�>�^>>�Tg>�~�>���>���>�׏    =���=�%�=�:�=�E=紽>H��>�L>�=�>�}T>�� >�_$>�j >���    =�I�=��j=�J=ˑ�=���>pF�>�.s>��>���?	��>��G>ϐ�>�Ib                                                            =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}�=�ξ>+�}=��=��>���?Q.?W�H?Y��??{�>���>	]�>{$    =��>l�>.��>S�>[c�>��>�9�>�b�>�0?5�?�b>�:>�>�    =��M=�^�>-[�>��>��>���? �?��?
�?"j">绅>�z}>��Q    =��u=��>-9�>	�%>��>���??}�?͟?#��>���>��>���    =ѹ�=��>-w�>"�r>)O�>�Y�>�^Z?��?��?�F>�K>�`;>���    =���=�4�>,�>$�>��>��I?N!?��?"��?%!L>��m>� >��Q    =�W�=��g>,ף>g>!g>�	?Z�?M?��?Xq>�`>>��[>�wi    =��m=���>)=�4�=φ�>���?5o?=v�?@e?2z>�v�>C>U�    =�	=�si>"r�=�Op=�~�>�Y�?(��?1�b?4��?+�;>��>S�>e\~    =��=���>!��=�t=��>�u!?-9�?5�"?8�d?-\Z>���>H��>[O�    =�vj=��>/�=Ւ�=�N�>��? ��?*D�?-_Z?'�>�	1>^!!>o��    =�Ӎ=�Z�>l�=���=�!>��?��?"NS?%��?#�=>��B>s�>�[�    =�|:=���>#��=�P�=���>���?Np?'�u?*��?&��>�B�>m�>}P�    =��R=�^�>$mf>:C>�>�k�?��?��?��?�Q>�4Y>�z�>���    =�p�=��>(=�>!�>'�1>���>�-�? ��?��?�->�tB>�d+>�
�    =�6=�=1>(�y><��>CҌ>�?>���>��F>���?	�m?A@>�2>�q"    =��;=�n~>(T�>!e�>(�>���>���? �&?�?�Z>�u>��o>�So    =��2=��P>'�>$E�>*��>��3>�)4>�<\?�S?m>�=N>��>�|�    =���=�}�>)�>,�v>3v�>��~>���>�>��3?�c>�HN>��>���    =���> L>+C>=��>Dλ>��>��h>�U$>�aQ?
��?�>ǈ>ț     =��>	4�>-�m>L��>T1M>��>��>��<>�5?�?j�>קD>�R�    =��>us>.�>S�`>[q�>��Z>�%a>�On>�$?1??�d>�E�>�I�    >Z�> >0A)>Ur_>\�>�J�>��0>���>��J?;�?z>��>��    >_�>+�>3m]>XuA>_�$>��>�.b>� >���?Q(?@>�>ݓ�    >
�?>�>9��>^��>f�>��!>��7>���>�6t?|�?��>�n�>ܚ�    >\�>&�>G`�>kd�>r��>��b>�Ȕ>��L>�'4?�?�S>�6�>ړv    >%>>4oD>U��>y�>�<�>�X>��>�Rg>�jx?9!?̵>�٬>�jU    >-l�>=={>^��>���>�kE>�[>��y>���>ɥ�?t�?+y>�f3>��    >/w�>?a�>`��>��~>�o�>�T�>��F>���>�n�?�??@>��>��    >0�p>@{�>a��>�=U>���>��I>�.P>���>�֞?��?��>���>֙    >4��>D�>f�>�;D>��&>���>�	|>ß�>�ky?�/?��>�'>��    >?�o>P<[>q��>���>�d�>��>�2�>�@g>��3?�?��>��B>��N    >Wc>h�5>�� >�&>��b>�#�>��\>�r9>څ�?N2?,�>΂�>��q    >u?�>���>�1O>�b2>�#>��>��>�Yf>���?˅>��>ǖ�>ȁ6    >��Q>�Ӵ>�fY>��2>��t>�b7>�J�>�̱>��?)=>�u>�>�أ    >��>�:>���>��1>���>�T0>ۦ�>��8>�jh?>>���>��>�7    >���>�#U>��H>���>�>�>�>���>���>�d�?*�>��>��>���    >�,A>�G=>��>���>�=�>�%>�G�>�
�>���?�>�j>���>��    >�9�>�D�>���>�k�>�D>�>�y�>�'}>��+?�>�*u>��>���    >��>���>�l�>�s>���>Ͷ^>�q�>�>���?��>�Sk>��>�޾    >��>���>�G�>���>�wJ>͠e>��J>��>���?Ԃ>�e�>�-�>��    >��>�|g>��M>�D�>��
>�P�>�N>�k>�d�?�0>�3�>��>��-    >��>�`5>���>� �>��X>�=�>���>���>��?`>��p>���>�    >��]>���>���>�>��X>�x'>��E>���>���?2l>��>���>���    >���>��6>�11>�B�>���>��t>��>�GL>��p?�>��>��'>��o    >�/>���>�C>�=P>��G>�ĉ>��>�V>��^?�>�4>�� >���    >�vd>��~>��M>���>�^l>�4>�п>�00>�K�?8�>�~�>�<>�h�    >�d>�]�>��>�h�>��>Ӧd>�u>�+F? ��?l�>��>���>���    >���>��L>��>���>�e�>�	/>��>��n? ��?W>�BR>�ݍ>�%�    >�n�>�˚>��>���>�D`>ϐ|>�N&>��>���?��>�u2>�\ >��y    >��Y>���>���>�0>���>̱>��>�>��?��>䞺>���>�(    >�Z.>�d>�^�>��l>���>��>���>�W>�q ??g>�$�>��>��p    >�ny>�g�>�>�>��t>�ik>¥X>�Dc>�U�>��Q?2�>�C>��E>�)
    >{�W>�=�>��[>�T�>��>��%>�1U>�e>��j? ��>��J>�U2>���    >]��>_��>~:>���>�P�>��>�Z>�P�>��*>�7_>ٲ?>�(>�b%    >Z4�>[C�>y˩>�g�>��>���>��>�Nd>��>���>�*�>��T>���    >_,�>aj>�o>�X�>��>���>��>�>�J�>�k=>��>�7�>���    >gm>j�v>��.>�;�>��.>���>�y	>���>�p�>�Ŝ>�	�>��>�m�    >s4�>xg~>���>�5�>��^>��>�Y(>�1>�>��>ܴ�>�d�>���    >���>��>�È>�4�>��>��>��6>���>��?4p>�ı>���>���                                                            4X�4X�3���;��b<�-O<��<�6<Q�=��=�ė=�w<��<��`=8�N<�J<��;z��<�1�<��c<x��;��:��&                                                                                                                                                                =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�    =}=�<�>,o�=��=���>�p�?SII?Y��?[Ɔ?@��>��>��>�                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    ?333?333?��?���?���?��?���?�`�?Z_?i�?d:�?tz?�?w�.?�3�?���?��H?�z�?���?�J�?���?��c?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2�?�2f?�,�?333