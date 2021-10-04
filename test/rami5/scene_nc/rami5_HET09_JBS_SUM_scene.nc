CDF       
      column        layer      >   layer_interface    ?   
wavelength              title         EForest description for the RAMI5 J�rvselja Birch Stand (Summer) scene      scene_id      HET09_JBS_SUM      scene         J�rvselja Birch Stand (Summer)     summary      jThis file contains a description of a RAMI5 scene ready to be read in by the SPARTACUS-Surface
radiative transfer model. In order that this model can treat tree trunks, we treat these forests as
vegetated urban areas, where the trunks are treated as buildings. The radiative effects of other
woody material (branches) has been added to the vegetation properties.     institution       2European Centre for Medium-Range Weather Forecasts     creator_name      Robin Hogan    creator_email         r.j.hogan@ecmwf.int          
wavelength                 	long_name         
Wavelength     units         nm     comment       UA dummy wavelength entry of zero is added at the end, for which all albedos are zero.         8  �   surface_type                	long_name         Surface type   
definition        ?0: Flat
1: Forest
2: Urban
3: Woody forest (or vegetated urban)         �   cos_solar_zenith_angle                  	long_name          Cosine of the solar zenith angle   units         1           �   top_flux_dn_sw                     	long_name         Top-of-canopy downwelling flux     units         W m-2         8  �   top_flux_dn_direct_sw                      	long_name         %Top-of-canopy downwelling direct flux      units         W m-2         8     nlayer                  	long_name         Number of active layers         P   height                     	long_name         Height of layer interfaces     units         m         �  T   building_fraction                      	long_name         Trunk cover fraction   units         1      comment       dThis variable is computed analytically from the number of trees falling withing the analysis domain.      �  P   veg_fraction                   	long_name         Vegetation fraction    units         1      comment       �This variable combines tree crowns and bare branches, and is computed from a rasterization of the
analysis domain in which circular trees are placed.         �  H   building_scale                     	long_name         Trunk effective diameter   units         m      comment       �This is the diameter of tree trunks in an equivalent forest in which all tree trunks have the same
diameter, and have the same areal fraction and surface area as the trees in the real forest.       �  @   	veg_scale                      	long_name         Vegetation horizontal scale    units         m      comment       �This is the variable S used in Eq. 20 of Hogan et al. (GMD 2018) relating the perimeter length of
the vegetation and the areal fraction.      �  8   veg_extinction                     	long_name         !Vegetation extinction coefficient      units         m-1    comment       �This variable is the wavelength-independent extinction coefficient of the vegetated fraction of
each layer, combining the contribution from foliage and tree branches (but not trunks) assuming the
surfaces of both to be randomly oriented.         �  0   veg_contact_fraction                   	long_name         Vegetation contact fraction    units         1      comment       �This is the fraction of the perimeter of the vegetation in each layer that is in contact with a
tree trunk, rather than clear air.        �  (   ground_sw_albedo                   	long_name         Ground albedo      units         1         8      wall_sw_albedo                        	long_name         Trunk albedo   units         1        �  X   roof_sw_albedo                        	long_name         Trunk narrowing albedo     units         1      comment      BThe area coverage of tree trunks reduces with height as the trunks get narrower. This reduction is
represented in a stepwise fashion, with a constant trunk diameter within individual layers. This
variable represents the albedo of the small horizontal surface exposed between two layers with
different trunk cover fraction.       �  *�   
veg_sw_ssa                        	long_name         #Vegetation single scattering albedo    units         1      comment       4Scattering by vegetation is assumed to be isotropic.     �  8x   foliage_extinction                     	long_name         Foliage extinction coefficient     units         m-1    comment      wThis variable is the wavelength-independent extinction coefficient of the vegetated fraction of
each layer, but containing only the contribution from foliage and not branches or trunks. This is
not used by SPARTACUS-Surface, but can be used (with foliage_sw_ssa) after running
SPARTACUS-Surface to work out the amount of radiation absorbed by foliage and branches separately.       �  F   foliage_sw_ssa                        	long_name          Foliage single scattering albedo   units         1      comment      This variable is like veg_sw_ssa but for the foliage only (no branches). It is not used by
SPARTACUS-Surface, but can be used (with foliage_extinction) after running SPARTACUS-Surface to
work out the amount of radiation absorbed by foliage and branches separately.     �  G C�y�C�9�D�D&P D*c3D1FfD<��DP` DXffD�P D�P E� E	p      �?   ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�   >�    ?   ?�  ?�  @   @   @@  @`  @�  @�  @�  @�  @�  @�  @�  @�  A   A  A  A  A   A(  A0  A8  A@  AH  AP  AX  A`  Ah  Ap  Ax  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  ;[�M;ND�;M�;L[;K�5;J��;I�;H�;Hg;G�;F;E��;D�;D	Q;C8X;Bj;>��;;0*;6U;/�};+F';'L;#a];X;S; �;��;	ZB;6�:���:�]:���:��:��0:Ē�:���:��@:�-�:Vw�:+��:��9�Y�9�89i�t9)&8��8��(8_Iv8!��7��17�L�7��'71��6�_56��6��5`��4n��2�{�0 �@$��d            7���<��<�5�=U=<�=�~=;��=4(�=��=}�=f�=�`==���>�<�>�h�>���>�X�>��@>��d>��+>���>~�m>�ml>�u)>��>��/>�8>��>��>§>�?
��?]=? �L?#v? m&??L?�-?=�?m�?:�?I�?�h>��>��J>��t>���>�̊>�,�>8��>*��>$��>H�>�n=#r�=*<��]<��D<kk�7���>a�M>Z�u>Y�E>Y��>Y7#>X��>X��>Xk�>XR>X*�>W�>W�!>W`�>WL3>WE�>W4�>U�^>T�R>SX4>Q�>O�M>O[c>Nh�>M�>M>M��>K,>I�>F�y>E�>D�A>E.C>F�,>F�>D6><��>2�6>&�2>��>
0�=���=��
=��[=���=��8=�|9=��*=|7H=�@=���=�}0=�Y-=t�=`Њ=:�7=��<�֢<,pJ;s�a:o�/�W\            ==�@#��@>��@Q`W@Y�?@[[H@f�%@O�@0�-@C3_@A�@Q%�@5�@sk�@|d�@|��@l��@h�@k�5@lj@_E@Te@>"=@I�-@K��@J
�@E��@Kg�@?�@IM�@(b�@+�,@0.�@5Xk@8T+@7h�@7��@6��@;�M@=�=@=�@>��@?t�@D|�@NF�@M�@KN�@G�!@A��@j�X@y�@vX�@ls�@i�K@��6@�l&@k�I@U߶@-�Y>���        >�{v<��k=��=��U=�=�|�=��>�u>�=��=�T=��E> /�=�3=���>;y�>_�>wR�>\��>e`�>g-�>{	p>u�0>ahM>��,>���>�;�>��$>|��>��>l�>l��>s#�>��:>��J>��>���>�� >�o>ŦB>�sK>���>�g >��>�H=>�DC>��K>H��>"��>�O >��9>n�>L��>0_�>���>�n�>���>a^�=Y�?<=�        =��_=+$A=,�=G<��f<��A<��=63=*�"<��b=)@�=4�=ȐC=%�=e|=l=>�=<��%<撊=v�D=oU�===m3<���<��F<��(=#��=��=��==2-=$&r=�=)$<�<��<�1#<´�<�x�<���<�O�<~�!<X�n<6=A<�;���;�͓;M��;*�2;�x;�(�;U�;$�:���;�F�;$ߠ:��9-n.3;ǰ    <s+�<��h=9�e<�o�<�~k=�x1>DN>U��>]�F>=�=��=*^x=DV    =Ǻ
=�>R�>��>��>tR�>��y>Û�>��y?o.>��>ϵ9>�C    =�=��>�>1v>>sѡ>�_�>�?6>�{6?i\>�:�>�=�>��p    =��|=�ט>�	>�>�U>s�<>�m:>�P(>Ό~?h�>�&�>�1V>��H    =�˿=Կ�>�t>�>�R>s�L>�s >�W�>Δ%?h~>��>�+�>��V    =ŵ�=ԣ[>��>�>�)>s�>�}�>�d�>Ρv?i8>��>�&�>�    =Ō#=�o�>�f>�@>.>s�.>��1>�v�>γ�?i>���>��>ϸa    =�i�=�D >�o>rC>Q>s� >���>Ĉ�>���?i�>��'>�>ϰh    =�A�=��>t�>C>�>s��>��.>Ě�>��L?i�>�܍>��>Ϧ�    =�-^=��>Z}>�>�!>s�^>���>ı�>��?l�>��,>�N>ϣ�    =��=��v>=G>� >�j>s��>��/>�� >�	?oh>���>��>ϟ�    =�*=ӷ�>&N>�4>��>sz.>���>��>�?r�>��>��>Ϟ�    =� =Ӭ>K>��>{�>sv�>���>��>�)�?t�>��8>�u>Ϟ7    =��=ӗ4>
�>�o>^�>sqd>��L>��5>�7m?v>���>� .>ϛ�    =��=ӆE>�0>j<>7�>slf>�3>�>�MV?y�>���>�.>Ϝ    =��V=�_>أ>8,>S>sb>� �>�'S>�e�?|g>��;>��u>ϗ    =ĕ=��>��>�c>�a>sO0>�<>�H�>χ�?~>�w�>��>ϊ�    =��=�F)>$\>+>�>s�>�t+>ŏ�>���?}5>�*�>�­>�a�    =�S=с�>�>`}>�>r�>���>���>�%2?��>��O>Ν�>�>�    =¾f=н>�>�->:�>r��>��M>�2�>�u ?>��U>�p�>��    =�\=�ף> �>��>HY>rgX>�7`>��>��N?u�>��>�0>�֍    =�2�=λ�=���>�H>0�>r�>��b>��h>�?nz>���>��>Γi    =���=�@<=��f>Fs>��>q�>��B>�W/>ѝ	?j4>��>͐0>�C�    =���=��=�
<>p>{e>qZ�>�:�>���>��?d�>�>�<+>���    =��G=�nR=�G!>��>$�>p��>���>�:�>҄�?_>���>��>ͨl    =�r�=��[=�j
>
b>�d>p��>��>ȵ5>�F?Yt>�jQ>̊>�U�    =��p=�$�=�9�>ɰ>�>p%e>�h�>�FT>ӕ?S�>�>� z>��$    =�K=Ĳ=�0s>��>
��>om�>��%>���>�Fh?9h>�>�p�>�T*    =���=��=��>��>	�>nϩ>�C�>�^@>ԯ?>ￆ>���>��%    =���=��%=�I�>�y>��>nG>�m�>ʕ]>���?��>���>�E�>�:    =��=�L/=��}>��>�>m��>��3>��O>�N?��>�*h>ɹ�>ʵT    =��H=��]=�ߛ>@*>7�>m-^>��>�>�S�?��>�P>�"+>�%�    =���=�C�=�כ=���>�~>l�\>��>�e�>ղs?n">�z�>ȓo>ɠ    =��=�E�=�XN=��>��>k��>�c�>��.>�@??�>�j�>�۳>��Z    =��s=�q{=��=��g> >k$�>���>�d>�e;? ��>�Gz>��>�2    =��@=�!�=�R�=�W=�|�>j^>���>�>>�X�? ��>� 0>�7 >�`5    =��=��p=��=��=�J2>i�C>�!�>ˊ~>��? f>��S>�9`>�ef    =�D=�\�=��=�8�=��>h��>�t�>ʶ<>��R>�ٯ>�w>�>�>G    =�v�=�"H=ᣵ=��=�>>g�>� F>�'
>�H>���>�r�>�=�>�l    =��R=��=�W�=�5�=���>g!�>�o�>�t�>ӊ#>�ɬ>�S^>�P�>À    =�L=���=��=���=�F�>fG�>���>ȟV>Ҧ�>�� >��>�EC>�t�    =��2=��9=��=�H�=��>e:�>�� >Ǫ�>ѡ�>���>ၸ>��u>�-�    =��)=��=�a=�c�=��\>c��>�Jy>��R>й
>�FW>߆>�j>��    =��=�/r=�E�=��_=��>b��>���>� e>���>���>ݲ�>���>�-�    =��X=�JV=��`=��=��>aq>��>�Ho>��>�#l>��<>���>�Î    =�9�=�.r=�qp=�=�Ԕ>`d�>��n>���>Θ0>��>�S�>�Q�>���    =��G=�=�/�=�S�=�w<>_��>�$�>�4�>��>��J>�G�>�t�>��z    =�^K=�Y$=�$�=�Ж=��@>^�t>��>���>�m>���>�#>�y@>���    =�(=�sr=�	2=��=��S>]�l>�Ӷ>3>�*I>�)�>�^Q>��>�^�    =���=�U�=ڐX=�{�=� >[3>�ZP>���>��>��>ҡ >��o>�*X    =�0~=�#v=�f�=���=��>V�D>���>��>�&>>�,A>�߻>�M$>���    =���=�U�=�=�=�Of=�r�>Q�>�f�>��_>��>�)>�a�>�8D>���    =��(=�e�=�Z�=��=�E>L��>��7>�30>���>���>�e�>�g�>���    =�5�=��=��=��=�G�>Go�>�j�>��%>��>�<�>�VX>���>�چ    =��=��&=��=�;=�.�>A1Y>�a�>�uy>�!�>��>�_>��>��@    =��Q=�t=�2�=�t�=�1>=�J>���>��(>�A>��>��]>�">�Q�    =��=���=�B�=�-=�K�>=eh>�x&>��h>��>���>�xq>���>��    =��=��=��u=�	�=���>>]9>�x�>�ܼ>�R�>�J�>��Q>�"�>�T�    =Ơa=��=�#=���=�u(>@�z>�Ĕ>��>�S�>�x>�/>��>�+e    =�܏=���=Ւ�=��=��>E�>�Z�>�^>>�Tg>�~�>���>���>�׏    =���=�%�=�:�=�E=紽>H��>�L>�=�>�}T>�� >�_$>�j >���    =�I�=��j=�J=ˑ�=���>pF�>�.s>��>���?	��>��G>ϐ�>�Ib                                                            =Ǻ
=�>R�>��>��>tR�>��y>Û�>��y?o.>��>ϵ9>�C    =�=��>�>1v>>sѡ>�_�>�?6>�{6?i\>�:�>�=�>��p    =��|=�ט>�	>�>�U>s�<>�m:>�P(>Ό~?h�>�&�>�1V>��H    =�˿=Կ�>�t>�>�R>s�L>�s >�W�>Δ%?h~>��>�+�>��V    =ŵ�=ԣ[>��>�>�)>s�>�}�>�d�>Ρv?i8>��>�&�>�    =Ō#=�o�>�f>�@>.>s�.>��1>�v�>γ�?i>���>��>ϸa    =�i�=�D >�o>rC>Q>s� >���>Ĉ�>���?i�>��'>�>ϰh    =�A�=��>t�>C>�>s��>��.>Ě�>��L?i�>�܍>��>Ϧ�    =�-^=��>Z}>�>�!>s�^>���>ı�>��?l�>��,>�N>ϣ�    =��=��v>=G>� >�j>s��>��/>�� >�	?oh>���>��>ϟ�    =�*=ӷ�>&N>�4>��>sz.>���>��>�?r�>��>��>Ϟ�    =� =Ӭ>K>��>{�>sv�>���>��>�)�?t�>��8>�u>Ϟ7    =��=ӗ4>
�>�o>^�>sqd>��L>��5>�7m?v>���>� .>ϛ�    =��=ӆE>�0>j<>7�>slf>�3>�>�MV?y�>���>�.>Ϝ    =��V=�_>أ>8,>S>sb>� �>�'S>�e�?|g>��;>��u>ϗ    =ĕ=��>��>�c>�a>sO0>�<>�H�>χ�?~>�w�>��>ϊ�    =��=�F)>$\>+>�>s�>�t+>ŏ�>���?}5>�*�>�­>�a�    =�S=с�>�>`}>�>r�>���>���>�%2?��>��O>Ν�>�>�    =¾f=н>�>�->:�>r��>��M>�2�>�u ?>��U>�p�>��    =�\=�ף> �>��>HY>rgX>�7`>��>��N?u�>��>�0>�֍    =�2�=λ�=���>�H>0�>r�>��b>��h>�?nz>���>��>Γi    =���=�@<=��f>Fs>��>q�>��B>�W/>ѝ	?j4>��>͐0>�C�    =���=��=�
<>p>{e>qZ�>�:�>���>��?d�>�>�<+>���    =��G=�nR=�G!>��>$�>p��>���>�:�>҄�?_>���>��>ͨl    =�r�=��[=�j
>
b>�d>p��>��>ȵ5>�F?Yt>�jQ>̊>�U�    =��p=�$�=�9�>ɰ>�>p%e>�h�>�FT>ӕ?S�>�>� z>��$    =�K=Ĳ=�0s>��>
��>om�>��%>���>�Fh?9h>�>�p�>�T*    =���=��=��>��>	�>nϩ>�C�>�^@>ԯ?>ￆ>���>��%    =���=��%=�I�>�y>��>nG>�m�>ʕ]>���?��>���>�E�>�:    =��=�L/=��}>��>�>m��>��3>��O>�N?��>�*h>ɹ�>ʵT    =��H=��]=�ߛ>@*>7�>m-^>��>�>�S�?��>�P>�"+>�%�    =���=�C�=�כ=���>�~>l�\>��>�e�>ղs?n">�z�>ȓo>ɠ    =��=�E�=�XN=��>��>k��>�c�>��.>�@??�>�j�>�۳>��Z    =��s=�q{=��=��g> >k$�>���>�d>�e;? ��>�Gz>��>�2    =��@=�!�=�R�=�W=�|�>j^>���>�>>�X�? ��>� 0>�7 >�`5    =��=��p=��=��=�J2>i�C>�!�>ˊ~>��? f>��S>�9`>�ef    =�D=�\�=��=�8�=��>h��>�t�>ʶ<>��R>�ٯ>�w>�>�>G    =�v�=�"H=ᣵ=��=�>>g�>� F>�'
>�H>���>�r�>�=�>�l    =��R=��=�W�=�5�=���>g!�>�o�>�t�>ӊ#>�ɬ>�S^>�P�>À    =�L=���=��=���=�F�>fG�>���>ȟV>Ҧ�>�� >��>�EC>�t�    =��2=��9=��=�H�=��>e:�>�� >Ǫ�>ѡ�>���>ၸ>��u>�-�    =��)=��=�a=�c�=��\>c��>�Jy>��R>й
>�FW>߆>�j>��    =��=�/r=�E�=��_=��>b��>���>� e>���>���>ݲ�>���>�-�    =��X=�JV=��`=��=��>aq>��>�Ho>��>�#l>��<>���>�Î    =�9�=�.r=�qp=�=�Ԕ>`d�>��n>���>Θ0>��>�S�>�Q�>���    =��G=�=�/�=�S�=�w<>_��>�$�>�4�>��>��J>�G�>�t�>��z    =�^K=�Y$=�$�=�Ж=��@>^�t>��>���>�m>���>�#>�y@>���    =�(=�sr=�	2=��=��S>]�l>�Ӷ>3>�*I>�)�>�^Q>��>�^�    =���=�U�=ڐX=�{�=� >[3>�ZP>���>��>��>ҡ >��o>�*X    =�0~=�#v=�f�=���=��>V�D>���>��>�&>>�,A>�߻>�M$>���    =���=�U�=�=�=�Of=�r�>Q�>�f�>��_>��>�)>�a�>�8D>���    =��(=�e�=�Z�=��=�E>L��>��7>�30>���>���>�e�>�g�>���    =�5�=��=��=��=�G�>Go�>�j�>��%>��>�<�>�VX>���>�چ    =��=��&=��=�;=�.�>A1Y>�a�>�uy>�!�>��>�_>��>��@    =��Q=�t=�2�=�t�=�1>=�J>���>��(>�A>��>��]>�">�Q�    =��=���=�B�=�-=�K�>=eh>�x&>��h>��>���>�xq>���>��    =��=��=��u=�	�=���>>]9>�x�>�ܼ>�R�>�J�>��Q>�"�>�T�    =Ơa=��=�#=���=�u(>@�z>�Ĕ>��>�S�>�x>�/>��>�+e    =�܏=���=Ւ�=��=��>E�>�Z�>�^>>�Tg>�~�>���>���>�׏    =���=�%�=�:�=�E=紽>H��>�L>�=�>�}T>�� >�_$>�j >���    =�I�=��j=�J=ˑ�=���>pF�>�.s>��>���?	��>��G>ϐ�>�Ib                                                                                                                                                                            =��j>rz>.�>S��>[o>���>�$>�NH>�?1+?��>�Fg>�I�    =�N�=���>' 3=�=�m�>�OM?"�}?,R�?/>�?2c4?j>>�Y�>�Zf    =��O=���>j�=�g�=���>�~F?9��?B��?E3;?Eb�?�x>��t>�q�    =��=� M>!4=�7�=�F">��I?6��?@=�?B�?F�m? �>�3�>���    =�i=��1>'B=��=��>���?A$�?J&?LY�?K>	?�Y>��>�BS    =r~q=��k>y=�~�=�tD>�7�?E�_?N��?P�n?P��?&mO>˓:>���    ={,=���> �=�ŝ=�a>�T ?D�J?L��?O#�?E_�?)�>���>��    =w�C=���>�=�qg=�׀>���?<#?D��?GD?@+?Vx>�3�>��Q    =}��=�:V>z�=�<=��j>�R�?7��?@��?C	?9m�>�h>~1�>�ޡ    =���=��d>�=ɩ!=�z�>�mJ?(d?1E�?40�?.<�>��>w��>��x    =��=��q>�=�9j=Ņ7>���?0R�?9��?<Yu?9�?�>�K>��n    =��m=���>P�=��
=�^k>�O�?9Y�?B��?E �?A��?lo>���>�ŗ    =�>�=��L>� =��=�<>��1?4��?>e�?@�?CTO?w>�X�>΂�    ={e�=�]v>&4=�cW=��}>��?@�U?Jm�?L��?Os?*4D>ֻ�>胗    =yD.=��>�2=���=��>�,�?@�4?Jud?L��?Q�?/�>ߏ>�H     =r��=��5>٩=��K=���>��?B��?L#&?NI�?P2_?)K�>���>�;�    =s�=�0t>|�=��	=�y�>���?B">?K��?M��?PB?)�5>���>��    =oc=��K>< =�_�=�[>���?D�p?NS�?Pn�?R�5?,�P>�u�>�2�    =rU=���>q�=�ď=�x�>�v?D&W?MƋ?O�7?S!0?/*�>ܸ�>�<    =x��=�6)>��=��&=�'>��{?B��?L��?N�?R��?0�C>�ش>�]5    =~=�K>`�=��8=�^�>��?Aog?Ke�?M�C?R�?18\>�N>��~    =��$=�|>e�=�i=��>�${?>�"?H��?K:?P�z?/>�Ȟ>��9    =��H=�Û>�=Ņv=śn>���?7��?A��?D&?K��?,�1>���>���    =�ؙ=��>!X=��=�)�>�ܳ?6H~?@bB?B��?J�?+m>�6>�b�    =�Z|=�a>+�m=�e�=�+�>�^�?:�W?DR�?F��?L͒?+�>��%>�X	    =� �=�6�>5_�=���=�L>��?9��?B�<?EK~?K\?*��>��O>��    =�}=�IX>=�=��u=�d�>��?7��?@�*?C�?Io�?)�>�d5>��3    =��=�~�>@��>V�>�H>���?5��?>�?AL�?H@?(C$>���>��    =�k�=��>@}�>X�>��>�u�?6??(?A��?HAy?(4�>քj>�z�    =���=�`�>A��> ��>c>�m?7��?@��?C;?IB?(��>���>�'s    =�y%=��X>D��>�>l5>���?:�\?C�G?E��?J��?)$�>Շl>��    =ϐI=��>U�>�>
��>�>(??1�?GJc?Ie'?L�?)�X>��>��%    =��z>*>mV�>�$>θ>�%?A�?Ia?J��?L̞?)�e>�K�>���    > �>��>�ھ>$�>(�z>�L�?Dc?J\�?L.F?L�w?);\>�o�>�H    >2>��>��>)��>-Ɗ>�ي?EG$?KJ�?MH?L�)?)!�>�e�>��    >��>�
>��>'1#>*�@>��G?F�?L��?NDL?MȐ?)�>�[>��;    >g}>{w>���>&�}>*3�>�Ʉ?Gi�?Mb??O{?NU
?)��>�F�>�̮    >8�>_=>��>&@
>)�>�]�?H+�?NT?O��?N̟?*7>�*>�7    >�> s>��q>%�3>)4�>��?H�?Ns7?P �?O?)�.>Ս�>�)|    =��4>:{>��@>"8�>%_>�L?I5k?OM�?P�?O��?*<�>��>�    =��>D�>�u|>!m>$�>ґ/?Ii?O�?Q3X?O�2?*D�>���>�l    =���>�K>�}d>��>!��>њ}?I�r?P�?Q�=?P-`?*n�>Ԥ?>娔    =�v>	l�>�8�>P.>�~>�F�?J�[?PŸ?Rn�?P��?*��>�f>喿    =�ٛ>O>�0�>�B>_�>Х?KLt?Q}<?S�?Q3�?+8>��>�6�    =�>�@>}F[>w�>�>�w?L�`?S�?T��?R�4?+�(>���>��    =�&>Ts>�"G>�l>>���?N�?T?U�g?S9b?,�~>��>���    =�D>R&>��
>7@>��>Ӵ�?O�%?U��?WA�?T��?-��>���>��    =�>c>���>a>v�>��?Ok�?UY&?V��?T$?-�>��>�M    =圑> �O>|��>rI>��>���?O�#?U�I?WP�?T�?-?m>��>��8    =�i�>m�>~�{>UG>��>��h?O �?U	�?V�J?S��?-/>�߶>��n    =��=�?O>{�G>�!>��>Б�?OH�?U:�?V�T?T?-~�>׍�>�Ԅ    =��%=�AF>t>��>
�>�uW?N��?U,(?V��?T"9?-�I>�p	>��G    =�X@=���>f��> Ç> ��>�u�?L�I?S��?U>?S�?,��>�;m>���    =�P�=��>U\a=��=�=�>���?L��?TE�?U�?T�?-��>�{|>蔦    =�lV=�I>R�=��=�Ծ>�H�?N�?V�?W�,?U�R?.�=>�.�>�0    =���=�]>RB�=�==ҝ�>��?P�K?X&8?Y�$?W�P?0W�>�6�>��    =�f�=���>TF=�[�=��>�)�?Sk?Z=�?[�8?YZ?1�@>�g�>샋    =��=�?>Z��=׃R=�3q>���?T�?[�?]Hd?Z��?2��>٦�>�ڋ    =��S=��z>]�
=�6�=���>���?[�?b5t?c�r?`&�?6ّ>�+M>�Z/    =�+9=_>�)�=�E=ۯS>��7?mƠ?pmP?q9�?k��??wK>��&>��                <��=۶=\�N=���=��=��=�F�=ڈ�<�=8wx=���=�T!=�q=�Z�>�+>7y!>N�`>7'q><z�>;�>GS�>6Y�>%0�>N=N>J�>;">9�b>1>=�>+��>.��>0�{>?�s>R�>n�->~��>��m>�S�>�l>���>�7M>�3�>�{�>��E>q=�>PWN> s,>�C>`�>ZY�>Au!>&��>�T>���>�ř>�JA>I��=M��<=�                                                                                                                                                                        =b�@=|aY>"��=���=�t>�b�?Wk6?^�/?`�o?Nb�?�t>|k>���    =N <=e!=>;=sM�=ry>��?Z�~?b��?dQ�?Y�?%M>�"�>��r    =C��=Y��>G9=i)=dK�>�d?\5?d��?f.3?^Rd?)?>�R>��    =I�5=`p�>�9=o1�=l�x>�z�?[A[?cw^?e�?['S?!n�>�fs>��    =A��=V��>r�=f��=a/�>���?\��?en?f��?_|=?,�>��>�O�    =^�0=w��>!2}=��=��m>��?Xm?_�?aMx?Pg<?�>��l>�=    =]�t=v��> ��=�yz=�'�>�v?X(�?_��?au�?P�?�>���>�Zf    =i�b=�F>%K�=�k=�qL>��~?VL!?]�?_I�?J�<>�w�>\'L>zH�    =q4�=��y>(�=�$=���>���?U!�?\/?]�K?F�a>�Q>:�T>U��    =\k=u�> C�=��T=�u>��?Xk?`4?a��?Q�R?
%>���>�h    =Y��=r|>{�=p=��>��{?X��?`�?b��?T�7?��>�T�>���    =P��=h��>m�=w�(=t��>��E?Zk�?c8?d��?[N?"��>�ͫ>��Z    =II=aa>6=q��=j(>��a?[A�?dv�?f?aY?2�G>�\,>� �    =B�=Y�>�:=km=aS>�ǃ?\&�?ec}?f��?c��?8�>��^>���    =A�=W
�>+�=g��=`C\>�[�?\n�?e#�?f�e?a?1 s>���>�$�    =@\�=V'>�s=f�/=_:h>�7�?\�/?eX�?f�?a`?1�>��g>�D;    =?�=US�>Y�=f)o=]��>��`?\׫?e��?gJ�?b�?4��>�5>���    =>f�=T�>��=eC�=[�>�`:?]7?f9?g�y?c�q?7`a>�՟>��    =@k�=V�t>m\=hY=^�>��c?\�a?f�?g��?d_s?9Z>�'O>�i    =B�=YFs>�=j�/=`K	>��1?\ע?f3W?g�5?d�?:��>�>�>��>    =@Q+=Vh�>��=g�-=]�H>��S?]?f3�?g�H?d��?9��>��x>�63    =@]=Uݛ>[�=g]�=]7�>��H?]I�?fl�?g��?d��?9op>� K>�k�    =Iz^=a�>�m=r�c=h��>�8�?\�&?e��?g��?dg�?8k5>��>�;    =X��=t(�>"��=��{=L+>��?]˺?f��?h14?d��?7�>�k�>���    =`tB=~��>)��=�=��>��u?^�>?g,#?h�q?e	�?7ͱ>ظ�>�*�    =k=�=�=>0�=�A�=��+>�HQ?_R�?g|+?i �?e?�?7��>�cg>�yp    =n�	=�GL>2��=�� =�u�>���?_P?gr�?h�r?e;�?7r`>���>�й    =n�=�?�>2=�=��[=�[ >��W?_C/?go\?h��?e7�?7Jy>�^�>�E�    =pg�=�,>2�G=�{�=�aK>��?_L-?gv?h�-?e;�?7/�>��>���    =q6�=��J>3��=�=��>�v?_q�?g�q?i�?eL?7E�>�$�>���    =r�=�c>CL�=��O=�s�>��a?a��?h��?jY�?fGj?8~A>���>�w    =�Y�=�^t>W'=��`=��]>�I�?ds?j��?k� ?gz-?9�J>�Ҁ>�0[    =���=��>i�=�p=�xD>���?f�v?l(�?mQr?h�f?;
>�Av>�L�    =�ݸ=�I'>p��=���=�
�>�i1?g�?l�d?m��?h��?;��>�o�>�b�    =��1=�a�>n$:=���=�:1>�T�?gR�?lH�?mio?h�]?;�$>߆`>���    =��%=�~�>o|z=�+�=�r�>��?geq?lNb?mm�?h��?;��>ߊ�>�x�    =�x�=��6>p�#=�ph=õ�>�Qi?gfT?lH�?mh�?h��?;��>�S�>�,�    =��>=��*>p�=���=è8>צ�?g.�?l"o?mEl?h��?;ix>��>���    =�1�=��>l��=��v=��>��A?f{b?k�\?l�V?h2�?:��>��e>�    =�k�=�L�>k��=�ד=�u>�Ir?fQ:?k� ?l��?h�?;)>���>��;    =�i=�3�>j�=���=�z>���?e��?kLT?l|T?g�s?:�b>���>�    =�w=��:>g�H=�=�S&>�!�?e�?j��?l*�?g�h?:�N>�v�>�F�    =���=���>h,�=��P=�R>�}�?e��?j�c?l,�?g��?:�>�">��    =��S=�?>eƭ=���=�A�>ϧ ?e�?j��?k��?gt!?:��>݃�>�\�    =�,^=�Y>h�~=��=�׆>��7?e�?k �?l/�?g��?;3�>��>���    =� #=�v�>k!�=�C=�q>Ӹ�?f�?kJ�?ls�?g�?;{>ߗn>��R    =��s=�"�>i��=��T=��/>��]?e�?k*C?lWo?gΈ?;;W>���>��?    =���=��u>g��=�r
=�P�>�5�?e|�?j�g?l0?g�8?;H>�s�>�X�    =�^�=�5>j�=��?=���>�F�?e�g?k�?lDh?g�/?;��>ࢷ>���    =���=��K>i��=�sa=��>�j�?e�d?j�?l$?g��?;��>�z>��    =�^=�[�>d�=��=��y>θ�?d��?j)�?kX�?g�?;��>�Et>�uK    =�_�=�ҙ>Ysc=���=��8>���?br�?h��?i�X?e�3?;K>�6}>�{    =y��=�{�>M�l=�.d=���>���?`�?fګ?h8;?d��?:16>�g�>��    =w��=�B�>K�j=��x=��<>�,?_��?f�t?g�;?dj�?:	�>�!>�o    =w��=�c�>K��=��=���>�V?_��?f�?g�\?dpn?:'>�!�>�w=    =y��=��>>M�V=�E)=��w>�ѽ?`�?f߆?h<�?d�)?:3�>�mB>��M    =�'�=�4w>S��=�%�=��>�f�?aF�?g�,?i�?e;�?:�|>�WT>��'    =��E=��>X>=��=��n>��C?b.4?hT�?i��?e�?;x>�	�>�T�    =�+9=_>�)�=�E=ۯS>��7?mƠ?pmP?q9�?k��??wK>��&>��    