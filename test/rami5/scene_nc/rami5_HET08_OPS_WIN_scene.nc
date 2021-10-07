CDF       
      column        layer      >   layer_interface    ?   
wavelength              title         CForest description for the RAMI5 Ofenpass Pine Stand (Winter) scene    scene_id      HET08_OPS_WIN      scene         Ofenpass Pine Stand (Winter)   summary      jThis file contains a description of a RAMI5 scene ready to be read in by the SPARTACUS-Surface
radiative transfer model. In order that this model can treat tree trunks, we treat these forests as
vegetated urban areas, where the trunks are treated as buildings. The radiative effects of other
woody material (branches) has been added to the vegetation properties.     institution       2European Centre for Medium-Range Weather Forecasts     creator_name      Robin Hogan    creator_email         r.j.hogan@ecmwf.int    software_version      0.7.3            
wavelength                 	long_name         
Wavelength     units         nm     comment       UA dummy wavelength entry of zero is added at the end, for which all albedos are zero.         8  �   surface_type                	long_name         Surface type   
definition        ?0: Flat
1: Forest
2: Urban
3: Woody forest (or vegetated urban)         �   cos_solar_zenith_angle                  	long_name          Cosine of the solar zenith angle   units         1           �   top_flux_dn_sw                     	long_name         Top-of-canopy downwelling flux     units         W m-2         8  �   top_flux_dn_direct_sw                      	long_name         %Top-of-canopy downwelling direct flux      units         W m-2         8  ,   nlayer                  	long_name         Number of active layers         d   height                     	long_name         Height of layer interfaces     units         m         �  h   building_fraction                      	long_name         Trunk cover fraction   units         1      comment       dThis variable is computed analytically from the number of trees falling withing the analysis domain.      �  d   veg_fraction                   	long_name         Vegetation fraction    units         1      comment       �This variable combines tree crowns and bare branches, and is computed from a rasterization of the
analysis domain in which circular trees are placed.         �  \   building_scale                     	long_name         Trunk effective diameter   units         m      comment       �This is the diameter of tree trunks in an equivalent forest in which all tree trunks have the same
diameter, and have the same areal fraction and surface area as the trees in the real forest.       �  T   	veg_scale                      	long_name         Vegetation horizontal scale    units         m      comment       �This is the variable S used in Eq. 20 of Hogan et al. (GMD 2018) relating the perimeter length of
the vegetation and the areal fraction.      �  L   veg_extinction                     	long_name         !Vegetation extinction coefficient      units         m-1    comment       �This variable is the wavelength-independent extinction coefficient of the vegetated fraction of
each layer, combining the contribution from foliage and tree branches (but not trunks) assuming the
surfaces of both to be randomly oriented.         �  D   veg_contact_fraction                   	long_name         Vegetation contact fraction    units         1      comment       |This is the fraction of the perimeter of the trunks in each layer that is in contact with
vegetation, rather than clear air.      �  <   ground_sw_albedo                   	long_name         Ground albedo      units         1         8  4   wall_sw_albedo                        	long_name         Trunk albedo   units         1        �  l   roof_sw_albedo                        	long_name         Trunk narrowing albedo     units         1      comment      BThe area coverage of tree trunks reduces with height as the trunks get narrower. This reduction is
represented in a stepwise fashion, with a constant trunk diameter within individual layers. This
variable represents the albedo of the small horizontal surface exposed between two layers with
different trunk cover fraction.       �  *�   
veg_sw_ssa                        	long_name         #Vegetation single scattering albedo    units         1      comment       4Scattering by vegetation is assumed to be isotropic.     �  8�   foliage_extinction                     	long_name         Foliage extinction coefficient     units         m-1    comment      wThis variable is the wavelength-independent extinction coefficient of the vegetated fraction of
each layer, but containing only the contribution from foliage and not branches or trunks. This is
not used by SPARTACUS-Surface, but can be used (with foliage_sw_ssa) after running
SPARTACUS-Surface to work out the amount of radiation absorbed by foliage and branches separately.       �  F   foliage_sw_ssa                        	long_name          Foliage single scattering albedo   units         1      comment      This variable is like veg_sw_ssa but for the foliage only (no branches). It is not used by
SPARTACUS-Surface, but can be used (with foliage_extinction) after running SPARTACUS-Surface to
work out the amount of radiation absorbed by foliage and branches separately.     �  GC�y�C�9�D�D&P D*c3D1FfD<��DP` DXffD�P D�P E� E	p      �?   ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�   >�    ?   ?�  ?�  @   @   @@  @`  @�  @�  @�  @�  @�  @�  @�  @�  A   A  A  A  A   A(  A0  A8  A@  AH  AP  AX  A`  Ah  Ap  Ax  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  A�  ;�l�;�bS;�آ;{�5;t\y;f;S��;Q7�;B�t;-�~;W;�m;��; �;
��:�(:�ݐ:�ϔ:�a�:|nr:Q��:`9��19_38�z.81+�7���77�5�H�                                                                                                                                    <�rs<���<P)<�%�=9��=G�P=v��=��=�Z�=�m=�ڟ=�j}=��>+�>- >&�>&M>/��>&�
>�#>�> ��=ׅk=���=d��=�<��]<9�>;ο�;"(                                                                                                                                >�co>�+�>��@>�Z>�V>��0>�V?>�' >}�>rHX>kZ�>kӨ>la�>k��>l�>]ă>RY�>L3U>EB4>BYu>K̚>$�=��8=� 1=��=�k0=��=Y�<��                                                                                                                                    ?tW�?v��?��?��:?�M_?�n-?�A�?�k?��E?��(?�b?��?��?��F?�EX?�C�?��?�(@�@$�@2*@�(@�@$�?�O�?�}'?�B�?��I?�(�?L�                                                                                                                                @R�#@�[?�K>���>��>dbH>qtF>�,�>�¡>��>��>��g>/�0>�%>	D�>e]�>r�>Q��>]1>}�>��>��n>���>���>��q?3�>�v>�s?G9|?�                                                                                                                                >k@= >�>WA�>��>��$>�z�>ڜx>ힽ?�h?

?��?Dc�?^��?b�8?g��?gd?f;?g�?ll?t��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�                                                                                                                                      ?a�a?e١?h1j?i�?i�<?g'$?f�?[MI?Y�>��=9�h<���=Z(5    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��    =�K=�53=���>!�N>+;2>M�z>uzx>�o�>��`?�]? �>�/@>��                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =��&=�Ւ>Tf�=���=���>�"�?R�?X�?Y�I?B\�>��j>#�|>:*�    =�,W=���>N�=��8=�|>���?J
�?P%�?R,�?>�>���>0�j>F��    =�l=���>Fnh=�-=���>��?=�e?D��?GB�?:
�>�&>C]�>Xt;    =���=�b<>C�*=��=�P�>�-C?9��?@�_?Cu�?8V�>ۀV>I��>^��    =��=�Ҍ>?�Z=�z.=�-]>���?3ޟ?;�"?>U�?6
�>�u�>R��>gA    =��_=��x>D}H=�ʈ=��o>�7?:��?B*:?D�?8�>��>G�P>\�    =�L�=�>HZm=�`�=�u>��-?@�R?G|�?I�?;/�>�>>�>TI�    =�W@=�}>H��=��6=Ϙ�>���?A#{?G�?J5w?;\�>��>>Q>S��    =�j�=�4>I3U=�1"=���>���?A��?H�r?J�t?;�>أ�>= >Rs<    =�˗=��r>K�=�a�=ʟ>�?E��?Lfa?N��?=O�>�A�>6�>L��    =��w=�T!>Jw=���=�r�>�T�?C<�?I֨?L?<47>�3�>;�>P��    =�	=���>?"�=�J=�>���?3$?:ʵ?=�U?5�u>��[>S��>hE    =�S�=�9>:'y=�D�=�Z>�?Q?+�Y?3��?6�?2��>�I�>_US>sQ    =�H=��>?!=�e�=�-%>�� ?2�?:�l?=v ?5��>��l>T)>ho�    =�$�=��>G:�=��=���>���??N?E��?HR�?:��>٤*>A��>V��    =���=��5>E�"=��I=���>�C?<��?C��?Fa�?9��>�b0>D�|>Y�9    =��]=�j@>C�9=���=��7>�ay?9�c?A25?C��?8wt>�dZ>Id|>^0j    =���=�E�>B��=�%_=�y�>�t ?8v�??�q?Bq�?7�\>��>K��>`O    =���=�X�>CPO=�p!=׵>���?96�?@��?C�?8/n>ۡ�>Jw�>_6�    =���=���>Eg=˂�=ԇ�>���?<QC?ClP?E��?9m�>ڑ�>E��>Z�    =��6=�^�>Cx�=�76=�wG>�Y?9s!?@â?CT ?8G�>ۍ(>JZ>^�v    =��=�=�>B�=�r2=��/>�@?8%??��?B)?7��>��t>L�>`��    =���=��X>D�:=�u*=Վ�>�F?;Pk?B~,?D� ?9>���>G?|>\%�    =�F(=��>H*�=ǣ�=�S�>��?@l�?G;)?I��?;�>�*=>?it>T��    =���=���>LD=�#A=�[>�Bs?FBB?L�?N�x?=jn>�+5>6w1>L,�    =��=�T�>J�=���=�l�>�X�?CB�?I�d?L�?<6�>�1�>;>P��    =���=��A>K�d=º]=��>��(?E� ?Ll?N;?=*->�b	>7l�>M�    =��&=��>S�n=��Y=�V>�O"?QM�?V�q?X��?A�>�c�>%��><�    =��=��>U=���=��N>��X?S��?Y?Z�7?B�>қD>"�>8�                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    @R�#@�C>��d>�\>]��><`�>M#�>�a>��->oD >���>��N>D=�O=�R�>C��>K&Y>+�F>30�>OU>r|�>y]�>���>��K>�GC?�>�D~>�|�?A��?/�                                                                                                                                =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\    =�9�=�L>V7�=��0=���>���?UL?Z��?\?�?Cr>��>g�>6:\                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    