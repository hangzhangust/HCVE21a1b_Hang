function a = cal_autocorrelation(E,step)
a = zeros(step,1);
    for i =1:step
        a(i) = (mean(E(:,1).*E(:,i+1),1)-mean(E(:,1))*mean(E(:,i+1)))/(mean(E(:,1).^2)-mean(E(:,1)).^2)^0.5/(mean(E(:,i+1).^2)-mean(E(:,i+1)).^2)^0.5;
    end
end