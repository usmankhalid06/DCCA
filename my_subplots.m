function my_subplots(nA,S,v,w,TCcorr,SMcorr,rTC,rSM)
    N = size(rTC{1},1);
    axis off
    set(gca,'Units','normalized','Position',[0 0 1 1]);
    str_title_h={'(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)'};
    str_title_v={'ST1','ST2','ST3','ST4','ST5','ST6','ST7','ST8'};

    

    a_vec = 0.96:-0.117:0.01;
    b_vec = [0.09 0.235 0.38 0.52 0.665 0.80 0.945]; %0.14
    c_vec = [0.165 0.305 0.445 0.59 0.73 0.87];
    d_vec = [0.225  0.365 0.51 0.65 0.79 0.93];
    for i =1:8
        text(0.005,a_vec(i), str_title_v{i},'Color','b','FontSize',12)
        if i<8
            text(b_vec(i), 0.015, str_title_h{i},'Color','b','FontSize',12)
        end
        if i>=2 && i<8
            text(c_vec(i-1),0.04, ['$\bar{\gamma}$ =' num2str(round(sum(SMcorr(i,:))/S,2),'%0.2f') ],'Color','k','FontSize',12, 'Interpreter','latex');
            text(d_vec(i-1),0.04, ['$\bar{\rho}$ =' num2str(round(sum(TCcorr(i,:))/S,2),'%0.2f') ],'Color','k','FontSize',12, 'Interpreter','latex')
        end

    end

    centerX = v/2; % X-coordinate of the circle center (adjust as needed)
    centerY = v/2; % Y-coordinate of the circle center (adjust as needed)
    radius = v/2;
    [x, y] = meshgrid(1:w, 1:v);
    distanceFromCenter = sqrt((x - centerX).^2 + (y - centerY).^2);
    circleMask = distanceFromCenter <= radius;
       

    ihs = 0.03;  %initial_horizontal_shift
    ivs = 1.00;  %initial_vertical_shift
    shz = 0.27/nA; %subplot_horizontal_size
    svs = S;  %subplot_vertical_size (more the better)
    vs  = 4.3/nA;  %horizontal shift of subplots
    nC  = S+0.55; %No. of rows
    shifter = 0.23;

    
    for i =1:nA 
        for j=1:S
            if i==1
            %%
            zscore_rxSM = abs(zscore(rSM{i}(j,:)));
            activationImage = flipdim(reshape(zscore_rxSM,v,w),3);
            maskedActivationImage = activationImage .* circleMask;
            maskedActivationImage(~circleMask) = 0;
            hax=axes();
            set(gca, 'Color', 'k'); hold on;
            imagesc(maskedActivationImage);
            cmap = jet(256);
            cmap(1,:) = 0;
            colormap(cmap);
            rectangle('Position', [centerX - radius, centerY - radius, 2 * radius, 2 * radius], ...
                'Curvature', [1, 1], 'EdgeColor', 'white', 'LineWidth', 1);

            newPos=[vs*(mod(j-1,1)+ihs+0.00),   (1-1/nC)-(1/nC)*(fix((j-1)/1)+ivs-1.2),   shz,   1/(svs+2)];
            set(gca,'outer',newPos),
            set(gca,'XTickLabel','')
            set(gca,'YTickLabel','')

            hax=axes();
            plot(zscore(rTC{i}(:,j)));  axis([0 N -3 3]);
            newPos=[vs*(mod(j-1,1)+ihs+0.04),   (1-1/nC)-(1/nC)*(fix((j-1)/1)+ivs-1.1),   3*shz,   1/(svs+1)];
            set(gca,'outer',newPos),

            else

            %%
            zscore_rxSM = abs(zscore(rSM{i}(j,:)));
            activationImage = flipdim(reshape(zscore_rxSM,v,w),3);
            maskedActivationImage = activationImage .* circleMask;
            maskedActivationImage(~circleMask) = 0;

            hax=axes();
            set(gca, 'Color', 'k'); hold on;
            imagesc(maskedActivationImage);
            cmap = jet(256);
            cmap(1,:) = 0;
            colormap(cmap);
            rectangle('Position', [centerX - radius, centerY - radius, 2 * radius, 2 * radius], ...
                'Curvature', [1, 1], 'EdgeColor', 'white', 'LineWidth', 1);
            newPos=[vs*(mod(j-1,1)+ihs+shifter*(i-1)),   (1-1/nC)-(1/nC)*(fix((j-1)/1)+ivs-1.2),   shz,   1/(svs+2)];
            set(gca,'outer',newPos),
            set(gca,'XTickLabel','')
            set(gca,'YTickLabel','')
            xlabel(['\gamma =' num2str(round(SMcorr(i,j),2))],'color','r','FontSize',11)
            xh = get(gca,'xlabel'); % handle to the label object
            p = get(xh,'position'); % get the current position property
            p(2) = 1+ p(2);       % double the distance,
            set(xh,'position',p)   % set the new position


            hax=axes();
            plot(zscore(rTC{i}(:,j))); axis([0 N -3 3]);
            newPos=[vs*(mod(j-1,1)+ihs+shifter*(i-1)+0.04),   (1-1/nC)-(1/nC)*(fix((j-1)/1)+ivs-1),   3*shz,   1/svs];
            set(gca,'outer',newPos),
            set(gca,'XTickLabel','')
            set(gca,'YTickLabel','')
            xlabel(['\rho',' = ',num2str(round(TCcorr(i,j),2))],'color','r','FontSize',11)
            end
        end
    end


