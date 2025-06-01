%% INTRODUCING STRAIGHT AND INCLINED CENTRAL CRACK
%% DATE: 31/07/2021

% % Intrducing a central crack of length 2ao
% Lexi = 0.5*M;
% Lexj = ceil(0.5*N)+1;
% Lex = (Lexi-1)*N + Lexj + 1;
% CrackNodes = Lex - ao + 1:Lex + ao - 1;
% for ID = Lex - ao + 1:Lex + ao - 1
%     for s = 1:length(node(ID).neighbors)
%         if(node(ID).aph(s) == 60)
%             node(ID).Kn(s) = 0;
%             node(ID).Kt(s) = 0;
%             node(ID).Kphi(s) = 0;
%             node(node(ID).neighbors(s)).Kn(node(node(ID).neighbors(s)).aph == -120) = 0;
%             node(node(ID).neighbors(s)).Kt(node(node(ID).neighbors(s)).aph == -120) = 0;
%             node(node(ID).neighbors(s)).Kphi(node(node(ID).neighbors(s)).aph == -120) = 0;
%             break;
%         end
%     end
%     for s = 1:length(node(ID).neighbors)
%         if(node(ID).aph(s) == 120)
%             node(ID).Kn(s) = 0;
%             node(ID).Kt(s) = 0;
%             node(ID).Kphi(s) = 0;
%             node(node(ID).neighbors(s)).Kn(node(node(ID).neighbors(s)).aph == -60) = 0;
%             node(node(ID).neighbors(s)).Kt(node(node(ID).neighbors(s)).aph == -60) = 0;
%             node(node(ID).neighbors(s)).Kphi(node(node(ID).neighbors(s)).aph == -60) = 0;
%             break;
%         end
%     end
% end
% if(ao ~= 0)
%     ID = Lex - ao;
%     for s = 1:length(node(ID).neighbors)
%         if(node(ID).aph(s) == 60)
%             node(ID).Kn(s) = 0;
%             node(ID).Kt(s) = 0;
%             node(ID).Kphi(s) = 0;
%             node(node(ID).neighbors(s)).Kn(node(node(ID).neighbors(s)).aph == -120) = 0;
%             node(node(ID).neighbors(s)).Kt(node(node(ID).neighbors(s)).aph == -120) = 0;
%             node(node(ID).neighbors(s)).Kphi(node(node(ID).neighbors(s)).aph == -120) = 0;
%             break;
%         end
%     end
% end

% Introducing inclined crack of length 2a

Lexi = 0.5*M;
Lexj = ceil(0.5*N)+1;
Lex = (Lexi-1)*N + Lexj + 1;

ID = Lex;
for i = 1:ao
    for s = 1:length(node(ID).neighbors)
        if(node(ID).aph(s) == 120)
            node(ID).Kn(s) = 0;
            node(ID).Kt(s) = 0;
            node(ID).Kphi(s) = 0;
            node(node(ID).neighbors(s)).Kn(node(node(ID).neighbors(s)).aph == -60) = 0;
            node(node(ID).neighbors(s)).Kt(node(node(ID).neighbors(s)).aph == -60) = 0;
            node(node(ID).neighbors(s)).Kphi(node(node(ID).neighbors(s)).aph == -60) = 0;
        end
        if(node(ID).aph(s) == 180)
            node(ID).Kn(s) = 0;
            node(ID).Kt(s) = 0;
            node(ID).Kphi(s) = 0;
            node(node(ID).neighbors(s)).Kn(node(node(ID).neighbors(s)).aph == 0) = 0;
            node(node(ID).neighbors(s)).Kt(node(node(ID).neighbors(s)).aph == 0) = 0;
            node(node(ID).neighbors(s)).Kphi(node(node(ID).neighbors(s)).aph == 0) = 0;
        end
    end
    ID = node(ID).neighbors(node(ID).aph == 60);
end

ID = node(Lex).neighbors(node(Lex).aph == -120);
for i = 1:ao-1
    for s = 1:length(node(ID).neighbors)
        if(node(ID).aph(s) == 120)
            node(ID).Kn(s) = 0;
            node(ID).Kt(s) = 0;
            node(ID).Kphi(s) = 0;
            node(node(ID).neighbors(s)).Kn(node(node(ID).neighbors(s)).aph == -60) = 0;
            node(node(ID).neighbors(s)).Kt(node(node(ID).neighbors(s)).aph == -60) = 0;
            node(node(ID).neighbors(s)).Kphi(node(node(ID).neighbors(s)).aph == -60) = 0;
        end
        if(i ~= ao-1)
            if(node(ID).aph(s) == 180)
                node(ID).Kn(s) = 0;
                node(ID).Kt(s) = 0;
                node(ID).Kphi(s) = 0;
                node(node(ID).neighbors(s)).Kn(node(node(ID).neighbors(s)).aph == 0) = 0;
                node(node(ID).neighbors(s)).Kt(node(node(ID).neighbors(s)).aph == 0) = 0;
                node(node(ID).neighbors(s)).Kphi(node(node(ID).neighbors(s)).aph == 0) = 0;
            end
        end
    end
    ID = node(ID).neighbors(node(ID).aph == -120);
end
