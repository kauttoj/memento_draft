function save_text(filename,data)

dlmwrite(filename,data,'delimiter','\t','precision',6);

end
