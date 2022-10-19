fprintf("\npostprocess\n");

if plot_options(2)
    figure
    hold on
    axis equal
    hh = trisurf(cnc', nodes(1,:), nodes(2,:), nodes(3,:), Jabs_real);
    set(hh,'FaceColor','flat');
    colormap('jet');
    colorbar;
    quiver3(midpoint(1,:), midpoint(2,:), midpoint(3,:), ...
            real(J(1,:)), real(J(2,:)), real(J(3,:)), ...
            'k', 'linewidth', 1.5);
    title('Volume current density - real part [A/m^2]')

    figure
    hold on
    axis equal
    hh = trisurf(cnc', nodes(1,:), nodes(2,:), nodes(3,:), Jabs_imag);
    set(hh,'FaceColor','flat');
    colormap('jet');
    colorbar;
    quiver3(midpoint(1,:), midpoint(2,:), midpoint(3,:), ...
            imag(J(1,:)), imag(J(2,:)), imag(J(3,:)), ...
            'k');
    title('Volume current density - imaginary part [A/m^2]')
endif
if plot_options(3)
    figure
    hold on
    axis equal
    hh = trisurf(cnc', nodes(1,:), nodes(2,:), nodes(3,:), norm([Jabs_real;Jabs_imag],'cols'));
    set(hh,'FaceColor','flat');
    colormap('jet');
    colorbar;
    title('Volume current density [A/m^2]')
endif
if plot_options(4)
    figure
    hold on
    axis equal
    hh = trisurf(cnc', nodes(1,:), nodes(2,:), nodes(3,:), losses);
    set(hh,'FaceColor','flat');
    colormap('jet');
    colorbar;
    title('Volumetric loss density [W/m^3]');
endif

if plot_options(5)
    if mod(plot_options(5),2) == 1
        figure
        plot(x_line(1,:),real(BtotLine(1,:)),'b-',...
            x_line(1,:),real(BtotLine(2,:)),'r-',...
            x_line(1,:),real(BtotLine(3,:)),'g-',...
            x_line(1,:),imag(BtotLine(1,:)),'b--',...
            x_line(1,:),imag(BtotLine(2,:)),'r--',...
            x_line(1,:),imag(BtotLine(3,:)),'g--')
        legend('Bx real','By real','Bz real','Bx imag','By imag','Bz imag')
        str_title = sprintf("B at %d Hz",freq);
        title(str_title)
    else
        figure
        plot(x_line(1,:),real(BtotLine(1,:)))
        title('Bx real')
        figure
        plot(x_line(1,:),real(BtotLine(2,:)))
        title('By real')
        figure
        plot(x_line(1,:),real(BtotLine(3,:)))
        title('Bz real')
        figure
        plot(x_line(1,:),imag(BtotLine(1,:)))
        title('Bx imag')
        figure
        plot(x_line(1,:),imag(BtotLine(2,:)))
        title('By imag')
        figure
        plot(x_line(1,:),imag(BtotLine(3,:)))
        title('Bz imag')
    endif


    for k = fields_id
        if plot_options(5) >= 3
        figure
            if exist('BtotLine_old')
                data = load(str_freq(k,:));
                if k/3 <= 1
                    plot(x_line(1,:),real(BtotLine_old(k,:)),'b-',x_line(1,:),real(BtotLine(k,:)),'g-',data(:,1),data(:,2),'r-')
                elseif k/3 > 1
                    plot(x_line(1,:),imag(BtotLine_old(k-3,:)),'b-',x_line(1,:),imag(BtotLine(k-3,:)),'g-',data(:,1),data(:,2),'r-')
                endif
                legend('PEEC std','PEEC G','benchmark')
                title(titles(k,:))
            else
                data = load(str_freq(k,:));
                if k/3 <= 1
                    plot(x_line(1,:),real(BtotLine(k,:)),'b-',data(:,1),data(:,2),'r-')
                elseif k/3 > 1
                    plot(x_line(1,:),imag(BtotLine(k-3,:)),'b-',data(:,1),data(:,2),'r-')
                endif
                legend('PEEC','benchmark')
                title(titles(k,:))
            endif
        endif
    end
endif  
  

if plot_options(6)
    figure
    hold on
    axis equal
    hh = trisurf(cnc', nodes(1,:), nodes(2,:), nodes(3,:), real(Bort(3,:)));
    set(hh,'FaceColor','flat');
    colormap('jet');
    colorbar;
    title('Bz field real on plate')
    figure
    hold on
    axis equal
    hh = trisurf(cnc', nodes(1,:), nodes(2,:), nodes(3,:), imag(Bort(3,:)));
    set(hh,'FaceColor','flat');
    colormap('jet');
    colorbar;
    title('Bz field imag on plate')
endif
fprintf("\n--------------------------------------------------\n")